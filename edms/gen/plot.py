### plot.py ###
# Author: Marc Zepeda
# Date: 2024-08-05

# Import packages
import os
import math
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
import numpy as np
from adjustText import adjust_text
from scipy.stats import ttest_ind
from ..gen import io
from ..gen import tidy as t

# Supporting methods
def re_un_cap(input_string: str):
    ''' 
    re_un_cap(): replace underscores with spaces and capitalizes each word for a given string
        
    Parameters:
    input_string (str): input string
    
    Dependencies:
    '''
    output_string = input_string.replace('_', ' ').title()
    return output_string

def round_up_pow_10(number):
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

def round_down_pow_10(number):
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

def log10(series):
    ''' 
    log10: returns log10 of maximum value from series or 0
    
    series: series, list, set, or array with values

    Dependencies: numpy
    '''
    return np.log10(np.maximum(series, 1))

def move_dist_legend(ax, legend_loc: str,legend_title_size: int,legend_size: int, legend_bbox_to_anchor: tuple, legend_ncol: tuple):
    ''' 
    move_dis_legend(): moves legend for distribution graphs.
    
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
    handles = old_legend.legendHandles
    labels = [t.get_text() for t in old_legend.get_texts()]
    title = old_legend.get_title().get_text()
    ax.legend(handles,labels,loc=legend_loc,bbox_to_anchor=legend_bbox_to_anchor,
              title=title,title_fontsize=legend_title_size,fontsize=legend_size,ncol=legend_ncol)

def extract_pivots(df: pd.DataFrame, x: str, y: str, vars='variable', vals='value'):
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
              title:str,title_size:int,title_weight:str,
              x_axis:str,x_axis_size:int,x_axis_weight:str,x_axis_scale:str,x_axis_dims:tuple,x_ticks_rot:int,xticks:list,
              y_axis:str,y_axis_size:int,y_axis_weight:str,y_axis_scale:str,y_axis_dims:tuple,y_ticks_rot:int,yticks:list,
              legend_title:str,legend_title_size:int,legend_size:int,legend_bbox_to_anchor:tuple,legend_loc:str,legend_items:tuple,legend_ncol:int,show:bool):
    ''' 
    formatter(): formats, displays, and saves plots.

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
    xticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    yticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    show (bool, optional): show plot (Default: True)
    
    Dependencies: os, matplotlib, seaborn, io, re_un_cap(), & round_up_pow_10()
    '''
    # Define plot types
    scats = ['scat', 'line', 'line_scat']
    cats = ['bar', 'box', 'violin', 'swarm', 'strip', 'point', 'count', 'bar_strip', 'box_strip', 'violin_strip','bar_swarm', 'box_swarm', 'violin_swarm']
    dists = ['hist', 'kde', 'hist_kde','rid']
    heats = ['ht']
        
    if typ not in heats:
        # Set title
        if title=='' and file is not None: title=re_un_cap(".".join(file.split(".")[:-1]))
        plt.title(title, fontsize=title_size, fontweight=title_weight)
        
        # Set x axis
        if x_axis=='': x_axis=re_un_cap(x)
        plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight)
        if x!='':
            if df[x].apply(lambda row: isinstance(row, (int, float))).all()==True: # Check that x column is numeric
                plt.xscale(x_axis_scale)
                if (x_axis_dims==(0,0))&(x_axis_scale=='log'): plt.xlim(round_down_pow_10(min(df[x])),round_up_pow_10(max(df[x])))
                elif x_axis_dims==(0,0): print('Default x axis dimensions.')
                else: plt.xlim(x_axis_dims[0],x_axis_dims[1])
        if xticks==[]: 
            if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(rotation=x_ticks_rot,ha='center')
            else: plt.xticks(rotation=x_ticks_rot,ha='right')
        else: 
            if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(ticks=xticks,labels=xticks,rotation=x_ticks_rot, ha='center')
            else: plt.xticks(ticks=xticks,labels=xticks,rotation=x_ticks_rot,ha='right')

        # Set y axis
        if y_axis=='': y_axis=re_un_cap(y)
        plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight)
        if y!='':
            if df[y].apply(lambda row: isinstance(row, (int, float))).all()==True: # Check that y column is numeric
                plt.yscale(y_axis_scale)
                if (y_axis_dims==(0,0))&(y_axis_scale=='log'): plt.ylim(round_down_pow_10(min(df[y])),round_up_pow_10(max(df[y])))
                elif y_axis_dims==(0,0): print('Default y axis dimensions.')
                else: plt.ylim(y_axis_dims[0],y_axis_dims[1])
        if yticks==[]: plt.yticks(rotation=y_ticks_rot)
        else: plt.yticks(ticks=yticks,labels=yticks,rotation=y_ticks_rot)

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

# Graph methods
def scat(typ: str,df: pd.DataFrame,x: str,y: str,cols=None,cols_ord=None,stys=None,cols_exclude=None,
         file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',
         figsize=(10,6),title='',title_size=18,title_weight='bold',
         x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,0),x_ticks_rot=0,xticks=[],
         y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,0),y_ticks_rot=0,yticks=[],
         legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),legend_ncol=1,show=True,
         **kwargs):
    ''' 
    scat(): creates scatter plot related graphs.

    Parameters:
    typ (str): plot type (scat, line, line_scat)
    df (dataframe): pandas dataframe
    x (str): x-axis column name
    y (str): y-axis column name
    cols (str, optional): color column name
    cols_ord (list, optional): color column values order
    stys (str, optional): styles column name
    cols_exclude (list, optional): color column values exclude
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    edgecol (str, optional): point edge color
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_scale (str, optional): x-axis scale linear, log, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    xticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    yticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    show (bool, optional): show plot (Default: True)
    
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), & round_up_pow_10()
    '''
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
              title,title_size,title_weight,
              x_axis,x_axis_size,x_axis_weight,x_axis_scale,x_axis_dims,x_ticks_rot,xticks,
              y_axis,y_axis_size,y_axis_weight,y_axis_scale,y_axis_dims,y_ticks_rot,yticks,
              legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol,show)

def cat(typ:str,df:pd.DataFrame,x='',y='',cols=None,cols_ord=None,cols_exclude=None,
        file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',lw=1,errorbar='sd',errwid=1,errcap=0.1,
        figsize=(10,6),title='',title_size=18,title_weight='bold',
        x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,0),x_ticks_rot=0,xticks=[],
        y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,0),y_ticks_rot=0,yticks=[],
        legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),legend_ncol=1,show=True, 
        **kwargs):
    ''' 
    cat(): creates category dependent graphs.

    Parameters:
    typ (str): plot type (bar, box, violin, swarm, strip, point, count, bar_swarm, box_swarm, violin_swarm)
    df (dataframe): pandas dataframe
    x (str, optional): x-axis column name
    y (str, optional): y-axis column name
    cols (str, optional): color column name
    cols_ord (list, optional): color column values order
    cols_exclude (list, optional): color column values exclude
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
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_scale (str, optional): x-axis scale linear, log, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    xticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    yticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    show (bool, optional): show plot (Default: True)
    
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), & round_up_pow_10()
    '''
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
        elif (x!='')&(y!=''): # x- and y-axis are specified
            if df[x].apply(lambda row: isinstance(row, str)).all()==True: # Check x colum is categorical
                cmap = cm.get_cmap(palette_or_cmap,len(df[x].value_counts()))
                palette = sns.color_palette([cmap(i) for i in range(cmap.N)])
            elif df[y].apply(lambda row: isinstance(row, str)).all()==True: # Check x colum is categorical
                cmap = cm.get_cmap(palette_or_cmap,len(df[y].value_counts()))
                palette = sns.color_palette([cmap(i) for i in range(cmap.N)])
        elif (x!='')&(y==''): # x-axis is specified
            if df[x].apply(lambda row: isinstance(row, str)).all()==True: # Check x colum is categorical
                cmap = cm.get_cmap(palette_or_cmap,len(df[x].value_counts()))
                palette = sns.color_palette([cmap(i) for i in range(cmap.N)])
                # Add palette = palette
        elif (x=='')&(y!=''): # y-axis is specified
            if df[y].apply(lambda row: isinstance(row, str)).all()==True: # Check x colum is categorical
                cmap = cm.get_cmap(palette_or_cmap,len(df[y].value_counts()))
                palette = sns.color_palette([cmap(i) for i in range(cmap.N)])
                # Add palette = palette
        else: return
    else: 
        print('Seaborn color palette or matplotlib color map not specified. Used seaborn colorblind.')
        palette = 'colorblind'

    fig, ax = plt.subplots(figsize=figsize)

    if cols is not None:

        if typ=='bar': sns.barplot(data=df, x=x, y=y, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='box': sns.boxplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='violin': sns.violinplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='swarm': sns.swarmplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='strip': sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='point': sns.pointplot(data=df, x=x, y=y, errorbar=errorbar, errwidth=errwid, capsize=errcap, hue=cols, hue_order=cols_ord, palette=palette, ax=ax, **kwargs)
        elif typ=='count': 
            if (x!='')&(y!=''):
                print('Cannot make countplot with both x and y specified.')
                return
            elif x!='': sns.countplot(data=df, x=x, hue=cols, hue_order=cols_ord, palette=palette, ax=ax, **kwargs)
            elif y!='': sns.countplot(data=df, y=y, hue=cols, hue_order=cols_ord, palette=palette, ax=ax, **kwargs)
            else:
                print('Cannot make countplot without x or y specified.')
                return
        elif typ=='bar_strip':
            sns.barplot(data=df, x=x, y=y, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='box_strip':
            sns.boxplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_strip':
            sns.violinplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='bar_swarm':
            sns.barplot(data=df, x=x, y=y, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='box_swarm':
            sns.boxplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_swarm':
            sns.violinplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        else:
            print('Invalid type! bar, box, violin, swarm, strip, point, count, bar_strip, box_strip, violin_strip, bar_swarm, box_swarm, violin_swarm')
            return

    else: # Cols was not specified
        
        if typ=='bar': sns.barplot(data=df, x=x, y=y, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='box': sns.boxplot(data=df, x=x, y=y, linewidth=lw, ax=ax, palette=palette, **kwargs)
        elif typ=='violin': sns.violinplot(data=df, x=x, y=y, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='swarm': sns.swarmplot(data=df, x=x, y=y, color=edgecol, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='strip': sns.stripplot(data=df, x=x, y=y, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='point': sns.pointplot(data=df, x=x, y=y, errorbar=errorbar, errwidth=errwid, capsize=errcap, palette=palette, ax=ax, **kwargs)
        elif typ=='count': 
            if (x!='')&(y!=''):
                print('Cannot make countplot with both x and y specified.')
                return
            elif x!='': sns.countplot(data=df, x=x, ax=ax, palette=palette, **kwargs)
            elif y!='': sns.countplot(data=df, y=y, ax=ax, palette=palette, **kwargs)
            else:
                print('Cannot make countplot without x or y specified.')
                return
        elif typ=='bar_strip':
            sns.barplot(data=df, x=x, y=y, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='box_strip':
            sns.boxplot(data=df, x=x, y=y, linewidth=lw, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_strip':
            sns.violinplot(data=df, x=x, y=y, edgecolor=edgecol, linewidth=lw, ax=ax, palette=palette, **kwargs)
            sns.stripplot(data=df, x=x, y=y, color=edgecol, edgecolor=edgecol, linewidth=lw, ax=ax, palette=palette, **kwargs)
        elif typ=='bar_swarm':
            sns.barplot(data=df, x=x, y=y, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='box_swarm':
            sns.boxplot(data=df, x=x, y=y, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_swarm':
            sns.violinplot(data=df, x=x, y=y, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        else:
            print('Invalid type! bar, box, violin, swarm, strip, point, count, bar_strip, box_strip, violin_strip, bar_swarm, box_swarm, violin_swarm')
            return

    formatter(typ,ax,df,x,y,cols,file,dir,palette_or_cmap,
              title,title_size,title_weight,
              x_axis,x_axis_size,x_axis_weight,x_axis_scale,x_axis_dims,x_ticks_rot,xticks,
              y_axis,y_axis_size,y_axis_weight,y_axis_scale,y_axis_dims,y_ticks_rot,yticks,
              legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol,show)

def dist(typ: str,df: pd.DataFrame,x: str,cols=None,cols_ord=None,cols_exclude=None,bins=40,log10_low=0,
        file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',lw=1,ht=1.5,asp=5,tp=.8,hs=0,des=False,
        figsize=(10,6),title='',title_size=18,title_weight='bold',
        x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,0),x_ticks_rot=0,xticks=[],
        y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,0),y_ticks_rot=0,yticks=[],
        legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),legend_ncol=1,show=True, 
        **kwargs):
    ''' 
    dist(): creates distribution graphs.

    Parameters:
    typ (str): plot type (hist, kde, hist_kde, rid)
    df (dataframe): pandas dataframe
    x (str): x-axis column name
    cols (str, optional): color column name
    cols_ord (list, optional): color column values order
    cols_exclude (list, optional): color column values exclude
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
    des (bool, optional): despine
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_scale (str, optional): x-axis scale linear, log, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    xticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    yticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    show (bool, optional): show plot (Default: True)
    
    Dependencies: os, matplotlib, seaborn, io, formatter(), re_un_cap(), & round_up_pow_10()
    
    Note: cannot set palette or cmap?
    '''
    # Omit excluded data
    if type(cols_exclude)==list: 
        for exclude in cols_exclude: df=df[df[cols]!=exclude]
    elif type(cols_exclude)==str: df=df[df[cols]!=cols_exclude]

    if typ=='hist':
        fig, ax = plt.subplots(figsize=figsize)
        if isinstance(bins, int):
            if x_axis_scale=='log': bins = np.logspace(log10(df[x]).min(), log10(df[x]).max(), bins + 1)
            else: bins = np.linspace(df[x].min(), df[x].max(), bins + 1)
        sns.histplot(data=df, x=x, kde=False, bins=bins, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        y=''
        y_axis='Count'
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        formatter(typ,ax,df,x,y,cols,file,dir,palette_or_cmap,
                  title,title_size,title_weight,
                  x_axis,x_axis_size,x_axis_weight,x_axis_scale,x_axis_dims,x_ticks_rot,xticks,
                  y_axis,y_axis_size,y_axis_weight,y_axis_scale,y_axis_dims,y_ticks_rot,yticks,
                  legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol,show)
    elif typ=='kde': 
        fig, ax = plt.subplots(figsize=figsize)
        if x_axis_scale=='log':
            df[f'log10({x})']=np.maximum(np.log10(df[x]),log10_low)
            sns.kdeplot(data=df, x=f'log10({x})', hue=cols, hue_order=cols_ord, linewidth=lw, ax=ax, **kwargs)
            x_axis_scale='linear'
            if x_axis=='': x_axis=f'log10({x})'
        else: sns.kdeplot(data=df, x=x, hue=cols, hue_order=cols_ord, linewidth=lw, ax=ax, **kwargs)
        y=''
        y_axis='Density'
        formatter(typ,ax,df,x,y,cols,file,dir,palette_or_cmap,
                  title,title_size,title_weight,
                  x_axis,x_axis_size,x_axis_weight,x_axis_scale,x_axis_dims,x_ticks_rot,xticks,
                  y_axis,y_axis_size,y_axis_weight,y_axis_scale,y_axis_dims,y_ticks_rot,yticks,
                  legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol,show)
    elif typ=='hist_kde':
        fig, ax = plt.subplots(figsize=figsize)
        if x_axis_scale=='log':
            df[f'log10({x})']=np.maximum(np.log10(df[x]),log10_low)
            bins = np.logspace(log10(df[x]).min(), log10(df[x]).max(), bins + 1)
            sns.histplot(data=df, x=f'log10({x})', kde=True, bins=bins, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
            x_axis_scale='linear'
            if x_axis=='': x_axis=f'log10({x})'
        else:
            bins = np.linspace(df[x].min(), df[x].max(), bins + 1) 
            sns.histplot(data=df, x=x, kde=True, bins=bins, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, ax=ax, **kwargs)
        y=''
        y_axis='Count'
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        formatter(typ,ax,df,x,y,cols,file,dir,palette_or_cmap,
                  title,title_size,title_weight,
                  x_axis,x_axis_size,x_axis_weight,x_axis_scale,x_axis_dims,x_ticks_rot,xticks,
                  y_axis,y_axis_size,y_axis_weight,y_axis_scale,y_axis_dims,y_ticks_rot,yticks,
                  legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol,show)
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
            ax.set_xlabel(x_axis,fontsize=x_axis_size,fontweight=x_axis_weight)
        g.set(yticks=yticks, ylabel=y_axis)
        g.set_titles("")
        if title=='' and file is not None: title=re_un_cap(".".join(file.split(".")[:-1]))
        g.figure.suptitle(title, fontsize=title_size, fontweight=title_weight)
        g.figure.subplots_adjust(top=tp,hspace=hs)
        if des==False: g.despine(top=False,right=False)
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

def heat(df: pd.DataFrame, x: str, y: str, vars='variable', vals='value',vals_dims:tuple=None,
         file=None,dir=None,edgecol='black',lw=1,annot=False,cmap="Reds",sq=True,cbar=True,
         title='',title_size=18,title_weight='bold',figsize=(10,6),
         x_axis='',x_axis_size=12,x_axis_weight='bold',x_ticks_rot=0,
         y_axis='',y_axis_size=12,y_axis_weight='bold',y_ticks_rot=0,
         show=True,**kwargs):
    '''
    heat(): creates heat plot related graphs.

    Parameters:
    df (dataframe): pandas dataframe
    x (str): x-axis column name
    y (str): y-axis column name
    vars (str, optional): variable column name
    vals (str, optional): value column name
    vals_dims (tuple, optional): vals minimum and maximum formatted (vmin, vmax; Default: None)
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
    figsize (tuple, optional): figure size per subplot
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_ticks_rot (int, optional): x-axis ticks rotation
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_ticks_rot (int, optional): y-axis ticks rotation
    show (bool, optional): show plot (Default: True)
    
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), & round_up_pow_10()
    '''
    # Find min and max values in the dataset for normalization
    if vals_dims is None:
        vmin = df[vals].values.min()
        vmax = df[vals].values.max()
    else:
        vmin = vals_dims[0]
        vmax = vals_dims[1]

    # Extract pivots
    dc = extract_pivots(df=df,x=x,y=y,vars=vars,vals=vals)

    # Create a single figure with multiple heatmap subplots
    fig, axes = plt.subplots(nrows=len(list(dc.keys())),ncols=1,figsize=(figsize[0],figsize[1]*len(list(dc.keys()))),sharex=False,sharey=True)
    if isinstance(axes, np.ndarray)==False: axes = np.array([axes]) # Make axes iterable if there is only 1 heatmap
    for (ax, key) in zip(axes, list(dc.keys())):
        sns.heatmap(dc[key],annot=annot,cmap=cmap,ax=ax,linecolor=edgecol,linewidths=lw,cbar=cbar,square=sq,vmin=vmin,vmax=vmax, **kwargs)
        if len(list(dc.keys()))>1: ax.set_title(key,fontsize=title_size,fontweight=title_weight)  # Add title to subplot
        else: ax.set_title(title,fontsize=title_size,fontweight=title_weight)
        if x_axis=='': ax.set_xlabel(re_un_cap(x),fontsize=x_axis_size,fontweight=x_axis_weight) # Add x axis label
        else: ax.set_xlabel(x_axis,fontsize=x_axis_size,fontweight=x_axis_weight)
        if y_axis=='': ax.set_ylabel(re_un_cap(y),fontsize=y_axis_size,fontweight=y_axis_weight) # Add y axis label
        else: ax.set_ylabel(y_axis,fontsize=y_axis_size,fontweight=y_axis_weight)
        plt.setp(ax.get_xticklabels(), rotation=x_ticks_rot, va='center', ha="right",rotation_mode="anchor") # Format x ticks
        plt.setp(ax.get_yticklabels(), rotation=y_ticks_rot, va='center', ha="right",rotation_mode="anchor") # Format y ticks
        ax.set_facecolor('white')  # Set background to transparent
    
    # Save & show fig
    if file is not None and dir is not None:
        io.mkdir(dir) # Make output directory if it does not exist
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight', format=f'{file.split(".")[-1]}')
    if show: plt.show()

def stack(df: pd.DataFrame,x:str,y:str,cols:str,cutoff=0,cols_ord=[],x_ord=[],
          file=None,dir=None,cmap='Set2',errcap=4,
          figsize=(10,6),title='',title_size=18,title_weight='bold',
          x_axis='',x_axis_size=12,x_axis_weight='bold',x_ticks_rot=45,x_ticks_ha='right',
          y_axis='',y_axis_size=12,y_axis_weight='bold',y_ticks_rot=0,
          legend_title='',legend_title_size=12,legend_size=12,
          legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_ncol=1,show=True,**kwargs):
    ''' 
    stack(): creates stacked bar plot

    Parameters:
    df (dataframe): pandas dataframe
    x (str, optional): x-axis column name
    y (str, optional): y-axis column name
    cols (str, optional): color column name
    cutoff (float, optional): y-axis values needs be greater than (e.g. 0)
    cols_ord (list, optional): color column values order
    cols_exclude (list, optional): color column values exclude
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    cmap (str, optional): matplotlib color map
    errcap (int, optional): error bar cap line width
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_ha (str, optional): x-axis ticks horizontal alignment
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_ticks_rot (int, optional): y-axis ticks rotation
    yticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    
    Dependencies: re, os, pandas, numpy, matplotlib.pyplot, & io
    '''
    # Make pivot table
    df_cut=df[df[y]>=cutoff]
    df_pivot=pd.pivot_table(df_cut, index=x, columns=cols, values=y, aggfunc=np.mean)
    df_pivot_err=pd.pivot_table(df_cut, index=x, columns=cols, values=y, aggfunc=np.std)
    if cols_ord!=[]: df_pivot=df_pivot[cols_ord]
    if x_ord!=[]: df_pivot=df_pivot.reindex(x_ord)

    # Make stacked barplot
    df_pivot.plot(kind='bar',yerr=df_pivot_err,capsize=errcap, figsize=figsize,colormap=cmap,stacked=True,**kwargs)

    # Set title
    if title=='' and file is not None: title=re_un_cap(".".join(file.split(".")[:-1]))
    plt.title(title, fontsize=title_size, fontweight=title_weight)
    
    # Set x axis
    if x_axis=='': x_axis=re_un_cap(x)
    plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight)
    plt.xticks(rotation=x_ticks_rot, ha=x_ticks_ha)
    
    # Set y axis
    if y_axis=='': y_axis=re_un_cap(y)
    plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight)
    plt.yticks(rotation=y_ticks_rot)

    # Set legend
    if legend_title=='': legend_title=re_un_cap(cols)
    plt.legend(title=legend_title, title_fontsize=legend_title_size, fontsize=legend_size, 
               bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol)
    
    # Save & show fig
    if file is not None and dir is not None:
        io.mkdir(dir) # Make output directory if it does not exist
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight', format=f'{file.split(".")[-1]}')
    if show: plt.show()

def cond_comp(df: pd.DataFrame, cond: str, cond_comp: str, wt:str, res: int, sample='sample', edit='edit', psuedocount=1):
    ''' 
    cond_comp(): returns DMS grid data in tidy format grouped by condition
    
    Need to update to be generalized to any dataset; test case Idris proteomics

    Parameters:
    df (dataframe): fastq outcomes dataframe
    cond (str): Condition column name for grouping fastq outcomes dataframe
    cond_comp (str): Condition for comparison group
    wt (str): Expected wildtype nucleotide sequence (in frame AA)
    res (int): First AA number
    sample (str, optional): Sample column name for splicing fastq outcomes dataframe (Default: 'sample')
    edit (str, optional): Edit column name within fastq outcomes dataframe (Default: 'edit')
    psuedocount (int, optional): psuedocount to avoid log(0) & /0 (Default: 1)
    
    Dependencies: Bio.Seq.Seq, pandas, numpy, tidy, edit_1(), dms_cond(), & aa_props
    '''
    #wt_nums = np.arange(res,res+len(wt_prot))
    print('Isolate single aa change fastq outcomes')
    dc=t.split((df),sample) # Isolate single aa change fastq outcomes and split by sample
    
    print('Fill with DMS grid data for each sample:')
    dc2=dict() # Fill with DMS grid data in tidy format split by sample
    for key_sample,df_sample in dc.items():
        print(key_sample)
        wt_fastq = df[(df['edit']=='WT')&(df[sample]==key_sample)] # Obtain WT fastq outcome
        df_sample_DMS=pd.DataFrame(columns=wt_fastq.columns) # Fill with DMS grid data in tidy format
        '''
        for num in wt_nums: # Iterate through WT protein sequence
            vals=dict() # Create dictionary with all amino acid changes for a given residue
            
            # Add metadata that is the same for all genotypes
            meta = [x for x in df_sample.columns if x not in [edit,'count','fraction','before','after','number']]
            for m in meta: 
                vals[m]=[wt_fastq[m].to_list()[0]]*len(list(aa_props.keys()))
            
            # Create all amino acid changes
            vals['before']=[wt_prot[num-res]]*len(list(aa_props.keys()))
            vals['number']=[num]*len(list(aa_props.keys()))
            vals['after']=list(aa_props.keys())
            vals[edit]=[vals['before'][i]+str(num)+vals['after'][i] for i in range(len(vals['after']))]

            # Fill in counts (+ psuedocount) for amino acid changes, WT, and none
            counts=[]
            num_mut = df_sample[df_sample['number']==num]
            for a in vals['after']:
                if a == wt_prot[num-res]: counts.append(wt_fastq['count'].to_list()[0]+psuedocount) # Wild type
                elif a in num_mut['after'].to_list(): counts.append(num_mut[num_mut['after']==a]['count'].to_list()[0]+psuedocount) # Amino acid change present
                else: counts.append(psuedocount) # Amino acid change absent
            vals['count']=counts
            sum_counts = sum(vals['count'])
            vals['fraction']=[count/sum_counts for count in vals['count']]

            df_sample_DMS = pd.concat([df_sample_DMS,pd.DataFrame(vals)]).reset_index(drop=True) # Append residue DMS data
        '''
        df_sample_DMS['number']=df_sample_DMS['number'].astype(int) # Set number as type int
        df_sample_DMS['count']=df_sample_DMS['count'].astype(int) # Set count as type int for plotting

        df_sample_DMS[sample] = [key_sample]*df_sample_DMS.shape[0]
        dc2[key_sample]=df_sample_DMS # Append sample DMS data

    print('Group samples by condition:')
    dc3=t.split(t.join(dc2,sample),cond) # Join samples back into 1 dataframe & split by condition
    df_cond_stat = pd.DataFrame()
    for key_cond,df_cond in dc3.items(): # Iterate through conditions
        print(key_cond)
        edit_ls = []
        fraction_avg_ls = []
        fraction_ls = []
        count_avg_ls = []
        before_ls = []
        after_ls = []
        number_ls = []
        for e in df_cond[edit]: # iterate through edits
            df_cond_edit = df_cond[df_cond[edit]==e]
            edit_ls.append(e)
            fraction_avg_ls.append(sum(df_cond_edit['fraction'])/len(df_cond_edit['fraction']))
            fraction_ls.append(df_cond_edit['fraction'].tolist())
            count_avg_ls.append(sum(df_cond_edit['count'])/len(df_cond_edit['count']))
            before_ls.append(df_cond_edit.iloc[0]['before'])
            after_ls.append(df_cond_edit.iloc[0]['after'])
            number_ls.append(df_cond_edit.iloc[0]['number'])
        df_cond_stat = pd.concat([df_cond_stat,
                                  pd.DataFrame({'edit':edit_ls,
                                                'before':before_ls,
                                                'after':after_ls,
                                                'number':number_ls,
                                                'fraction_ls':fraction_ls,
                                                'fraction_avg':fraction_avg_ls,
                                                'count_avg':count_avg_ls,
                                                cond:[key_cond]*len(number_ls)})])
    df_cond_stat = df_cond_stat.drop_duplicates(subset=['edit','Description']).reset_index(drop=True)

    # Fold change & p-value relative comparison group
    print(f'Compute FC & pval relative to {cond_comp}:')
    df_stat = pd.DataFrame()
    df_comp = df_cond_stat[df_cond_stat[cond]==cond_comp] # Isolate comparison group
    df_other = df_cond_stat[df_cond_stat[cond]!=cond_comp] # From other groups
    for e in set(df_other[edit].tolist()): # iterate through edits
        print(f'{e}')
        df_other_edit = df_other[df_other[edit]==e]
        df_comp_edit = df_comp[df_comp[edit]==e]
        df_other_edit['fraction_avg_compare'] = [df_comp_edit.iloc[0]['fraction_avg']]*df_other_edit.shape[0]
        df_other_edit['count_avg_compare'] = [df_comp_edit.iloc[0]['count_avg']]*df_other_edit.shape[0]
        df_other_edit['FC'] = df_other_edit['fraction_avg']/df_comp_edit.iloc[0]['fraction_avg']
        ttests = [ttest_ind(other_fraction_ls,df_comp_edit.iloc[0]['fraction_ls']) 
                                 for other_fraction_ls in df_other_edit['fraction_ls']]
        df_other_edit['pval'] = [ttest[1] for ttest in ttests]
        df_other_edit['tstat'] = [ttest[0] for ttest in ttests]
        df_stat = pd.concat([df_stat,df_other_edit])
    df_stat['compare'] = [cond_comp]*df_stat.shape[0]
    return df_stat[[edit,'before','after','number','FC','pval','tstat','fraction_avg','fraction_avg_compare','count_avg','count_avg_compare',cond,'compare']].sort_values(by=['number','after']).reset_index(drop=True) 

def vol(df: pd.DataFrame,x: str,y: str,stys:str=None, size:str=None,size_dims:tuple=None,
        file=None,dir=None,palette_or_cmap='YlOrRd',edgecol='black',
        figsize=(10,6),title='',title_size=18,title_weight='bold',
        x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_dims=(0,0),x_ticks_rot=0,xticks=[],
        y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_dims=(0,0),y_ticks_rot=0,yticks=[],
        legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',
        legend_items=(0,0),legend_ncol=1,display_size=True,display_labels=True,return_df=True,show=True,
        **kwargs):
    ''' 
    vol(): creates volcano plot
    
    Parameters:
    df (dataframe): pandas dataframe
    x (str): x-axis column name
    y (str): y-axis column name
    cols (str, optional): color column name
    sty (str, optional): sty column name
    size (str, optional): size column name
    size_dims (tuple, optional): (minimum,maximum) values in size column (Default: None)
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    edgecol (str, optional): point edge color
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    xticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    yticks (list, optional): y-axis tick values
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
    
    Dependencies: os, matplotlib, seaborn, pandas, & edit_1()
    '''
    # Strings with subscripts
    log2 = 'log\u2082'
    log10 = 'log\u2081\u2080'
    
    # Log transform data
    df[f'{log2}({x})'] = [np.log10(xval)/np.log10(2) for xval in df[x]]
    df[f'-{log10}({y})'] = [-np.log10(yval) for yval in df[y]]
    
    # Organize data by significance
    signif = []
    for (log2FC,log10P) in zip(df[f'{log2}({x})'],df[f'-{log10}({y})']):
        if (np.abs(log2FC)>1)&(log10P>-np.log10(0.05)): signif.append(f'FC & p-value')
        elif (np.abs(log2FC)<=1)&(log10P>-np.log10(0.05)): signif.append('p-value')
        elif (np.abs(log2FC)>1)&(log10P<=-np.log10(0.05)): signif.append('FC')
        else: signif.append('NS')
    df['Significance']=signif
    signif_order = ['NS','FC','p-value','FC & p-value']

    # Organize data by abundance
    sizes=(1,100)
    if size_dims is not None: df = df[(df[size]>=size_dims[0])&(df[size]<=size_dims[1])]

    # Set dimensions
    if x_axis_dims==(0,0): x_axis_dims=(min(df[f'{log2}({x})']),max(df[f'{log2}({x})']))
    if y_axis_dims==(0,0): y_axis_dims=(0,max(df[f'-{log10}({y})']))

    # Generate figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # with significance boundraries
    plt.vlines(x=-1, ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    plt.vlines(x=1, ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    plt.hlines(y=-np.log10(0.05), xmin=x_axis_dims[0], xmax=x_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    
    # with data
    if display_size==False: size=None
    sns.scatterplot(data=df, x=f'{log2}({x})', y=f'-{log10}({y})', 
                    hue='Significance', hue_order=signif_order, 
                    edgecolor=edgecol, palette=palette_or_cmap, style=stys,
                    size=size, sizes=sizes,
                    ax=ax, **kwargs)
    
    # with labels
    if display_labels:
        df_signif = df[df['Significance']=='FC & p-value']
        adjust_text([plt.text(x=df_signif.iloc[i][f'{log2}({x})'], 
                              y=df_signif.iloc[i][f'-{log10}({y})'],
                              s=edit) for i,edit in enumerate(df_signif['edit'])])

    # Set title
    if title=='' and file is not None: title=p.re_un_cap(file[-4])
    plt.title(title, fontsize=title_size, fontweight=title_weight)
    
    # Set x axis
    if x_axis=='': x_axis=f'{log2}({x})'
    plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight)
    if xticks==[]: 
        if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(rotation=x_ticks_rot,ha='center')
        else: plt.xticks(rotation=x_ticks_rot,ha='right')
    else: 
        if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(ticks=xticks,labels=xticks,rotation=x_ticks_rot, ha='center')
        else: plt.xticks(ticks=xticks,labels=xticks,rotation=x_ticks_rot,ha='right')

    # Set y axis
    if y_axis=='': y_axis=f'-{log10}({y})'
    plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight)

    if yticks==[]: plt.yticks(rotation=y_ticks_rot)
    else: plt.yticks(ticks=yticks,labels=yticks,rotation=y_ticks_rot)

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