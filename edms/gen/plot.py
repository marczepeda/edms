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
from ..gen import io

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
        if title=='' and file is not None: title=re_un_cap(file[-4])
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
        if title=='' and file is not None: title=re_un_cap(file[-4])
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

def heat(typ: str,df: pd.DataFrame, x: str, y: str, vars='variable', vals='value',
         file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',lw=1,annot=False,cmap="Reds",sq=True,
         title='',title_size=18,title_weight='bold',
         x_axis='',x_axis_size=12,x_axis_weight='bold',x_ticks_rot=0,
         y_axis='',y_axis_size=12,y_axis_weight='bold',y_ticks_rot=0,
         show=True,**kwargs):
    '''
    heat(): creates heat plot related graphs.

    Updates: 
    Need to make split heat maps optional
    Need to keep palette_or_cmap or cmap
    Need to add figsize maybe.
    Needs to be updated to be similar to edms.fastq.dms_grid()

    Parameters:
    typ (str): plot type (ht)
    df (dataframe): pandas dataframe
    x (str): x-axis column name
    y (str): y-axis column name
    vars (str, optional): variable column name
    vals (str, optional): value column name
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    edgecol (str, optional): point edge color
    lw (int, optional): line width
    annot (bool, optional): annotate values
    cmap (str, optional): matplotlib color map
    sq (bool, optional): square
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
    
    Note: cannot set palette or cmap?
    '''
    # Find min and max values in the dataset for normalization
    vmin = df[vals].values.min()
    vmax = df[vals].values.max()

    # Generates dictionary of pivot tables
    pivots = extract_pivots(df=df,x=x,y=y,vars=vars,vals=vals)
    pivots['colorbar']=pivots[next(iter(pivots))]

    # Create a single figure with multiple heatmap subplots
    fig, axes = plt.subplots(nrows=1,ncols=len(list(pivots.keys())),figsize=(pivots[next(iter(pivots))].shape[0]*len(list(pivots.keys())),pivots[next(iter(pivots))].shape[1]),sharex=True,sharey=True)
    for i, (ax, subtitle) in enumerate(zip(axes, list(pivots.keys()))):
        if i!=len(axes)-1: sns.heatmap(pivots[subtitle],annot=annot,cmap=cmap,square=sq,ax=ax,vmin=vmin,vmax=vmax,linecolor=edgecol,linewidths=lw,cbar=False) # No heatmap colorbar
        else: # Add heatmap colorbar to the end
            sns.heatmap(pivots[subtitle],annot=annot,cmap=cmap,ax=ax,vmin=vmin,vmax=vmax,linecolor=edgecol,linewidths=lw,cbar=True)
            plt.gca().set_visible(False)
        if len(list(pivots.keys()))>1: ax.set_title(subtitle,fontsize=x_axis_size,fontweight=x_axis_weight)  # Add title to subplot
        else: ax.set_title(title,fontsize=x_axis_size,fontweight=x_axis_weight)
        if i==0: # First heatmap subplot has x & y axes
            if x_axis=='': ax.set_xlabel(re_un_cap(x),fontsize=x_axis_size,fontweight=x_axis_weight)
            else: ax.set_xlabel(x_axis,fontsize=x_axis_size,fontweight=x_axis_weight)
            if y_axis=='': ax.set_ylabel(re_un_cap(y),fontsize=y_axis_size,fontweight=y_axis_weight)
            else: ax.set_ylabel(y_axis,fontsize=y_axis_size,fontweight=y_axis_weight)
            ax.tick_params(axis='x',rotation=x_ticks_rot)
            ax.tick_params(axis='y',rotation=y_ticks_rot)
        else: # Remove redundant information from remaining heatmap subplots 
            if x_axis=='': ax.set_xlabel(re_un_cap(x),fontsize=x_axis_size,fontweight=x_axis_weight)
            else: ax.set_xlabel(x_axis,fontsize=x_axis_size,fontweight=x_axis_weight)
            ax.set_ylabel('')
            ax.tick_params(axis='x',rotation=x_ticks_rot)
        ax.set_facecolor('white')  # Set background to transparent

    formatter(typ,ax,df,x,y,'cols',file,dir,palette_or_cmap,
              title,title_size,title_weight,
              x_axis,x_axis_size,x_axis_weight,'linear',(0,0),x_ticks_rot,
              y_axis,y_axis_size,y_axis_weight,'linear',(0,0),y_ticks_rot,
              '',12,9,(1,1),'upper left',(0,0),show)

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
    if title=='' and file is not None: title=re_un_cap(file[-4])
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

def vol(df: pd.DataFrame,x: str,y: str, stys=None,size=None,cols_exclude=None,
        file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',
        figsize=(10,6),title='',title_size=18,title_weight='bold',
        x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,0),x_ticks_rot=0,xticks=[],
        y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,0),y_ticks_rot=0,yticks=[],
        legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),legend_ncol=1,show=True,
        **kwargs):
    ''' 
    vol(): creates volcano plot
    
    Work in progress...

    Parameters:
    df (dataframe): pandas dataframe
    x (str): x-axis column name
    y (str): y-axis column name
    cols (str, optional): color column name
    stys (str, optional): styles column name
    size (str, optional): size column name
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
    
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), & round_up_pow_10()
    '''
    # Strings with subscripts
    log2 = 'log\u2082'
    log10 = 'log\u2081\u2080'
    
    # Log transform data
    df[f'{log2}({x})'] = [-np.log10(xval)/np.log10(2) for xval in df[x]]
    df[f'-{log10}({y})'] = [-np.log10(yval) for yval in df[y]]
    
    # Organize data by significance
    signif = []
    for (log2FC,log10P) in zip(df[f'{log2}({x})'],df[f'-{log10}({y})']):
        if (np.abs(log2FC)>1)&(log10P>-np.log10(0.05)): signif.append(f'{log2}FC & p-value')
        elif (np.abs(log2FC)<=1)&(log10P>-np.log10(0.05)): signif.append('p-value')
        elif (np.abs(log2FC)>1)&(log10P<=-np.log10(0.05)): signif.append(f'{log2}FC')
        else: signif.append('NS')
    df['Significance']=signif
    signif_order = ['NS',f'{log2}FC','p-value',f'{log2}FC & p-value']

    # Organize data by conservation

    # Set dimensions
    if x_axis_dims==(0,0): x_axis_dims=(min(df[f'{log2}({x})']),max(df[f'{log2}({x})']))
    if y_axis_dims==(0,0): y_axis_dims=(0,max(df[f'-{log10}({y})']))

    # Generate figure
    fig, ax = plt.subplots(figsize=(5,5))
    
    # with significance boundraries
    plt.vlines(x=-1, ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    plt.vlines(x=1, ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    plt.hlines(y=-np.log10(0.05), xmin=x_axis_dims[0], xmax=x_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    
    # and with data
    sns.scatterplot(data=df, x=f'{log2}({x})', y=f'-{log10}({y})', hue='Significance', hue_order=signif_order, edgecolor=edgecol, palette='YlOrRd', ax=ax, **kwargs)

    # Move legend to the right of the graph
    ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
              bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol)

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