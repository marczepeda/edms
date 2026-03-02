# COMMON IMPORTS #
import sys
import os
import gc
import math
import time
import random
import pandas as pd
import numpy as np
from pathlib import Path
import seaborn as sns
from biopandas.pdb import PandasPdb

# MATPLOTLIB #
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from matplotlib import gridspec

# DATA CLASS IMPORTS #
from dataclasses import dataclass, field, replace
from typing import Any, Dict, List, Optional, Tuple, Union
from pathlib import Path

# IGNORE WARNING #
import logging
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

# SCIPY AND STATS #
import scipy.cluster as sp_cl
import scipy.interpolate as interp
from scipy.stats import pearsonr
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster

import statsmodels.stats.multitest as smm
from statsmodels.nonparametric.smoothers_lowess import lowess
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

import pickle

# SETUP MPL SETTINGS #

mpl.rcParams.update({
    # Keep text as text in SVG
    'svg.fonttype': 'none',

    # Font
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',

    # Font sizes
    'font.size': 6,
    'axes.labelsize': 6,
    'axes.titlesize': 6,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,

    # Thin axes lines
    'axes.linewidth': 0.5,
    'lines.linewidth': 0.5,

    # Thin tick lines
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 2,
    'ytick.major.size': 2,

    # Smaller padding
    'axes.labelpad' : 2,
    'ytick.major.pad' : 1,
    'xtick.major.pad': 1
})

cm_2_in = 0.393
in_2_cm = 2.54
five_across_size = (4.3/in_2_cm, 2.8/in_2_cm)
three_across_size = (6.35/in_2_cm, 2.8/in_2_cm)

# Dataclass for axes, labels, etc
@dataclass(frozen = True)
class AxisLabelOpts:
    title: Optional[str] = None
    xlabel: Optional[str] = None
    ylabel: Optional[str] = None
    xwindow: Optional[Tuple[float, float]] = None
    xlim: Optional[Tuple[float, float]] = None
    ylim: Optional[Tuple[float, float]] = None
    xticks: Optional[List[int]] = None
    yticks: Optional[List[int]] = None
    protein_len: Optional[int] = None
    linewidth: Optional[float] = 0.5
    tick_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"labelsize":8, "width":0.5})
    space: Optional[float] = 3

@dataclass(frozen = True)
class NegCtrlOpts:
    ###
    adjust: bool = False
    neg_ctrl_col: Optional[str] = 'gene'
    neg_ctrl_conditions: Optional[str] = 'NT Controls'
    lines: bool = False
    sd: Optional[float] = 1
    ctrl_line_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"color": "k", "ls": "--", "lw": 1})

# Dataclass for style, colors, etc
@dataclass(frozen = True)
class ScatterStyleOpts:
    palette: Union[str, List[str]] = "colorblind"
    hue_col: Optional[str] = None
    marker_col: Optional[str] = None
    marker_map: Optional[Dict[str, str]] = field(default_factory=dict)
    color_map: Optional[Dict[Any, str]] = None
    transparent_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"edgecolor": "black", "alpha": 0.6, "s": 3, "linewidth": 0.1})
    opaque_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"edgecolor": "black", "alpha": 1, "s": 5, "linewidth": 0.25})
    highlight_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"edgecolor": "black", "alpha": 1, "s": 5, "linewidth": 0.3})
    transparency_threshold: Optional[float] = 2.0
    line_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"color": "k", "ls": "--", "lw": 1})
    cutoff_line_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"color": "gray", "ls": ":", "lw": 0.5})

@dataclass(frozen = True)
class BoxStyleOpts:
    color_map: Optional[Dict[Any, str]] = None
    xlabel_rotation: Optional[float] = 45
    xlabel_ha: Optional[str] = "right"

    axhline_kws: Dict[str, Any] = field(default_factory=lambda:
        {'linestyle': '--', 'color': 'grey'})
    box_kws: Dict[str, Any] = field(default_factory=lambda:
        {"width":0.4, "showfliers":False, "saturation":1})
    boxprops: Dict[str, Any] = field(default_factory=lambda:
        {'linewidth':0.5, 'edgecolor':'black'})
    whiskerprops: Dict[str, Any] = field(default_factory=lambda:
        {'linewidth':0.5})
    capprops: Dict[str, Any] = field(default_factory=lambda:
        {'linewidth':0.5})
    medianprops: Dict[str, Any] = field(default_factory=lambda:
        {'linewidth':1.5, 'color':'black'})

@dataclass(frozen = True)
class JitterStyleOpts:
    linewidth: Optional[float] = 0.5
    colors: Optional[List[str]] = field(default_factory=lambda:
        ['red', 'purple', 'blue', 'cornflowerblue'])
    axhline_kws: Dict[str, Any] = field(default_factory=lambda:
        {'linestyle': '--', 'color': 'grey'})

    box_kws: Dict[str, Any] = field(default_factory=lambda:
        {'width':0.3, 'showcaps':True, 'showfliers':False})
    stripplot_kws: Dict[str, Any] = field(default_factory=lambda:
        {'alpha': 0.75, 'size': 2, 'linewidth': 0, 'zorder':0,
         'jitter':True, 'dodge':False, 'rasterized':True})
    stripplot_hits_kws: Dict[str, Any] = field(default_factory=lambda:
        {'alpha': 1.0, 'size': 3, 'linewidth': 0.25, 'zorder':0,
         'jitter':True, 'dodge':False, 'rasterized':False})
    kdeplot_kws: Dict[str, Any] = field(default_factory=lambda:
        {'linewidth':0.5})

    boxprops: Dict[str, Any] = field(default_factory=lambda:
        {'facecolor': 'none', 'zorder': 1, 'linewidth': 0.5})
    whiskerprops: Dict[str, Any] = field(default_factory=lambda:
        {'color': 'gray', 'zorder': 1, 'linewidth': 0.5})
    capprops: Dict[str, Any] = field(default_factory=lambda:
        {'zorder': 1, 'color': 'gray', 'linewidth': 0.5})
    medianprops: Dict[str, Any] = field(default_factory=lambda:
        {'zorder': 1, 'color': 'gray', 'linewidth': 0.5})

@dataclass(frozen = True)
class LollipopStyleOpts:
    line_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"color":'black', "linewidth":0.5, "linestyle":'--'})
    marker_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"markeredgecolor":"black", "markersize":3, "markeredgewidth":0.4})
    stem_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"color":"black", "linewidth":0.4})
    base_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"color":"black"})
    marker_rasterized: bool = True
    stem_rasterized: bool = True

@dataclass(frozen = True)
class LoessStyleOpts:
    line_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"color": "k", "ls": "--", "lw": 1})
    pval_plot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"alpha": 1, "markersize": 1})
    regr_plot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"alpha": 1, "markersize": 1, "color":"gray"})
    below_color: Optional[str] = "lightgray"
    above_color: Optional[str] = "red"

@dataclass(frozen = True)
class ClusteredHeatmapStyleOpts:
    col_width: float = 0.05
    heatmap_width: float = 0.75
    wspace: float = 0.02
    color_pal: str = "tab20"
    default_color: str = "#EEEEEE"
    palette: Dict[str, Any] = field(default_factory=lambda: None)
    heatmap_kws: Dict[str, Any] = field(default_factory=lambda:
        {"cmap":"bwr", "center":0, "vmin":-10, "vmax":15, "cbar":False, "rasterized":True})
    border_kws: Dict[str, Any] = field(default_factory=lambda:
        {"fill":False, "clip_on":False, "zorder":10, "edgecolor":"black", "linewidth":0.75})
    line_kws: Dict[str, Any] = field(default_factory=lambda:
        {"color":'black', "linewidth":0.75})

@dataclass(frozen = True)
class PWESHeatmapStyleOpts:
    # HEATMAPS #
    heatmap_kws: Dict[str, Any] = field(default_factory=lambda:
        {"square":True, "cmap":'RdBu_r', "vmin":-1, "vmax":1,
         "xticklabels":False, "yticklabels":False, "cbar":False})
    clustered_heatmap_kws: Dict[str, Any] = field(default_factory=lambda:
        {"vmin":-1, "vmax":1, "center":0, "cmap":'RdBu_r',
         "xticklabels":False, "yticklabels":False})
    line_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"color": "k", "ls": "--", "lw": 0.25})

    swarmplot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"size": 2, "linewidth": 0, "edgecolor":None})
    jitterplot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"alpha":0.5, "size":2, "linewidth":0, "edgecolor":None})
    stripplot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"alpha":0.5, "size":2, "linewidth":0, "edgecolor":None})

    boxplot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"fliersize":0, "linewidth":0.5, "width":0.5, })
    boxprops: Dict[str, Any] = field(default_factory=lambda:
        {'facecolor': 'none', 'edgecolor': "black"})
    whiskerprops: Dict[str, Any] = field(default_factory=lambda:
        {'color': 'black'})
    capprops: Dict[str, Any] = field(default_factory=lambda:
        {'color': 'black'})
    medianprops: Dict[str, Any] = field(default_factory=lambda:
        {'color': 'black'})

    spine_color: str = "black"
    spine_visible: bool = True
    spine_wspace: float = 0.2
    spine_hspace: float = 0.2
    linewidth: float = 0.5
    edgecolor: str = "black"
    # COLOR BAR #
    colorbar_normalize_kws: Dict[str, Any] = field(default_factory=lambda:
        {"vmin":-1, "vmax":1})
    colorbar_base_kws: Dict[str, Any] = field(default_factory=lambda:
        {"cmap":'RdBu_r'})
    colorbar_linewidth: float = 0.5

@dataclass(frozen = True)
class DomainOpts:
    dom_setting: Dict[str, List[Dict[str, Any]]] = field(default_factory=dict)
    gene: str = field(default_factory=str)
    default_alpha: float = 0.25
    default_color: str = "lightblue"
    gene_col: str = "Gene Symbol"

@dataclass(frozen = True)
class LegendOpts:
    loc: str = None
    ncol: int = 1
    frameon: bool = False
    mark_size: Optional[int] = 5
    title: Optional[str] = None
    path: Union[str, Path] = "legend"

# Dataclass for outputs
@dataclass(frozen = True)
class OutputOpts:
    figsize: Tuple[float, float] = (4.0, 2.2)
    subplots_kws: Dict[str, Any] = field(default_factory=dict)
    save: bool = False
    path: str = "scatterplot"
    show: bool = True
    dpi: int = 1200
    transparent: bool = True
    rasterize: bool = True
    out_type: str = "svg"

# Defining default options
AXIS = AxisLabelOpts()
SCATTERSTYLE = ScatterStyleOpts()
NEGCONTROL = NegCtrlOpts()
LOLLIPOPSTYLE = LollipopStyleOpts()
BOXPLOTSTYLE = BoxStyleOpts()
JITTERPLOTSTYLE = JitterStyleOpts()
LOESSSTYLE = LoessStyleOpts()
CLUSTEREDHEATMAPSTYLE = ClusteredHeatmapStyleOpts()
PWESHEATMAPSTYLE = PWESHeatmapStyleOpts()
DOMAIN = DomainOpts()
LEGEND = LegendOpts()
OUTPUT = OutputOpts()

arial_font6 = fm.FontProperties(fname=f'Arial.ttf', size=6)
 
