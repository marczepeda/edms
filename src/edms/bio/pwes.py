'''
Module: fastq.py
Author: Marc Zepeda
Created: 2026-03-01
Adapted from: BE-SCAN (Calvin Hu)
Description: Proximity-Weighted Enrichment Score (PWES) clustering

Usage:
[Dataclasses] ### From BE-SCAN ###
- [Dataclass for axes, labels, etc]
    - AxisLabelOpts: Options for axis labels, ticks, limits, etc.
    - PWESHeatmapStyleOpts: Options for PWES heatmap styling.
- [Dataclass for outputs]
    - OutputOpts: Options for figure output (size, path, format, etc.)
    - Defining default options: Instances of the above dataclasses with default values.
[Helper Functions]
- [HEX <-> RGB CONVERSION]
    - hex_to_rgb: Convert hex color code to RGB tuple.
    - hex_list_to_rgb_list: Convert list of hex color codes to list of RGB tuples.
    - tuple_to_hex: Convert RGB tuple to hex color code.
    - fillgaps: Fill in missing amino acid positions in the PWES matrix with NaN values.
- [DISTANCE FACTOR (e^(-d^2 / 2t^2))]
    - process_pdb: Process PDB file to extract centroid coordinates for amino acids.
    - get_pairwise_dist: Calculate pairwise distances between amino acids from centroid coordinates.
    - gauss: Apply Gaussian distance function to pairwise distances.
- [ENRICHMENT SCORE CALCULATIONS]
    - hill: Apply Hill function to scale log2 fold change scores.
    - calculate_pw_score: Calculate pairwise sums matrix for scores and apply tanh transformation.
- [Calculate PWES]
    - calculate_pwes: Calculate PWES by multiplying pairwise score matrix with Gaussian distance matrix, and optionally zeroing out negative scores.
    - cluster_pwes: Perform hierarchical clustering on PWES matrix and assign cluster labels.
    - get_clus_aa: Print amino acids in each cluster from clustering results.
    - shuffle_pwes: Randomize PWES scores by shuffling and multiplying with Gaussian distance matrix, to generate a null distribution for comparison.
- [Iterative Clustering]
    - safe_slug: Make a filesystem-friendly folder name component.
    - mask_equals: Equality mask that treats NaN == NaN for grouping purposes.
### Test ###
- [PWES Cluster Mapping]
    - pymol_script(): Generates a PyMOL script to display selected residues as separate objects.

[Plot Functions]
- hist: Generates a histogram of the number edits for each cluster.
- cat: creates categorical graphs for PWES 3D clustering results
- torn: generates tornado plots for visualizing PWES clusters
- heatmap: generates heatmap of PWES scores, with options for clustering and showing colorbar
- heatmap_cbar: generates a standalone colorbar for the heatmap with specified orientation and styling.
### Test ###
- clustermap: generates clustered heatmap of PWES scores, with options for coloring clusters and styling.

[Main Function]
### Revise to include pymol_script() automation ### Test ###
- clustering: Compute 3D PWES clustering analysis
    - run_one: Inner function that runs the clustering pipeline for a single group of data (used in iterative mode)
'''
# Import packages
import sys
import os
import gc
import math
import time
import random
from importlib.resources import files
from datetime import datetime
import mpld3
import re
import pandas as pd
import numpy as np
from pathlib import Path
import seaborn as sns
from biopandas.pdb import PandasPdb
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.font_manager as fm
import matplotlib.patches as patches

from dataclasses import dataclass, field, replace
from pathlib import Path

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

from ..gen import io, tidy as t, plot as p, stat as st
from ..data import uniprot, pdb
from .fastq import add_label_info
from ..utils import mkdir

# Helper functions
### HEX <-> RGB CONVERSION
def hex_to_rgb(hex_color):
    """
    Converts a hex color code to an RGB tuple.
    """
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i+2], 16)/256 for i in (0, 2, 4))

def hex_list_to_rgb_list(hex_color_list):
    """
    Converts a list of hex color codes to a list of RGB tuples.
    """
    return [hex_to_rgb(hex_color) for hex_color in hex_color_list]

def tuple_to_hex(rgb_tuple):
    '''
    Converts an RGB tuple to a hex color code.
    '''
    if not (isinstance(rgb_tuple, tuple) and len(rgb_tuple) == 3 and
            all(isinstance(val, float) and 0 <= val <= 1 for val in rgb_tuple)):
        return None

    r, g, b = rgb_tuple
    r, g, b = int(r * 255), int(g * 255), int(b * 255)
    return f"{r:02x}{g:02x}{b:02x}".upper()

### FILL GAPS
def fillgaps(
    df_pwes_sorted, gene_bounds,
    ):

    labels_list = []
    for gene, (start, end) in gene_bounds.items():
        labels_list += [gene + str(i).zfill(4) for i in range(start, end + 1)]

    # Create DataFrame with NaNs
    full_matrix = pd.DataFrame(np.nan, index=labels_list, columns=labels_list)

    common_rows = full_matrix.index.intersection(df_pwes_sorted.index)
    common_cols = full_matrix.columns.intersection(df_pwes_sorted.columns)
    full_matrix.loc[common_rows, common_cols] = df_pwes_sorted.loc[common_rows, common_cols]

    return full_matrix

### DISTANCE FACTOR (e^(-d^2 / 2t^2))

def process_pdb(pdb_file, chains):
    '''
    Given PDB or Alphafold .pdb file, returns dataframe with X, Y, Z coordinates and amino acid number.

    pdb_file: str, path to pdb file
    '''
    ppdb = PandasPdb().read_pdb(pdb_file)
    atom_df = ppdb.df['ATOM'] # return just atoms
    if len(chains) > 0:
        atom_df = atom_df.loc[atom_df["chain_id"].isin(chains)]

    # gather XYZ of alpha carbons
    coord_df = atom_df.loc[atom_df["atom_name"] == "CA"]
    df_centroids = coord_df[["chain_id", "residue_number", "x_coord", "y_coord", "z_coord"]]
    df_centroids.columns = ["chain", "aa_num", "x", "y", "z"]

    if len(chains) > 0: df_centroids['label'] = df_centroids['chain'] + df_centroids['aa_num'].astype(str).str.zfill(4)
    else: df_centroids['label'] =  df_centroids['aa_num'].astype(str)

    return df_centroids

def get_pairwise_dist(df_centroids, chains, aa_int=None):
    """
    Calculate pairwise distances from centroid coordinates.
    Returns a pandas df of pairwise distances with columns/index as AA pos

    df_centroids: pandas dataframe containing 4 columns of
        ['aa_num', 'x', 'y', 'z'].
    aa_int: optional tuple of (aa_min, aa_max) defining the aas to calculate pdists
        the default is None, which takes the min/max of df_centroids
    """

    df_centroids['aa_num'] = df_centroids['aa_num'].astype('int64')

    # Isolate desired amino acid interval
    if aa_int is None:
        # remove unresolved residues (xyz = NaN) before finding aa min/max
        print("Removing unresolved residues...")
        df_aaint = df_centroids.loc[~df_centroids.isnull().any(axis=1)].copy()

    else:
        aa_min = aa_int[0]
        aa_max = aa_int[1]
        print("Removing unresolved residues...")
        df_aaint = df_centroids.loc[df_centroids['aa_num'].between(aa_min, aa_max)].copy()
        df_aaint = df_aaint.loc[~df_aaint.isnull().any(axis=1)].copy()
        if df_aaint['aa_num'].min() != aa_min:
            print('Warning! User aa_min input was ' + str(aa_min))
            print('But first resolved AA was ' + str(df_aaint['aa_num'].min()))
        if df_aaint['aa_num'].max() != aa_max:
            print('Warning! User aa_max input was ' + str(aa_max))
            print('But last resolved AA was ' + str(df_aaint['aa_num'].max()))

    # calculate all pairwise distances in euclidean 3d space, condense to square-form
    print("Calculating pairwise distances...")
    pairwise = pdist(df_aaint[['x','y','z']], 'euclidean')
    pairwise = squareform(pairwise)

    if len(chains) > 0: df_pwdist = pd.DataFrame(pairwise, index=df_aaint['label'], columns=df_aaint['label'])
    else: df_pwdist = pd.DataFrame(pairwise, index=df_aaint['aa_num'], columns=df_aaint['aa_num'])
    return df_pwdist

def gauss(distance, std):
    """
    Applies Gaussian distance function.
    """
    arg = -(distance * distance) / (2 * std * std)
    dist = np.exp(arg)
    return dist

### ENRICHMENT SCORE CALCULATIONS

def hill(lfc, m, theta):
    """
    From NZL code:
    hill function to scale the log2_fc numbers into range [0,1]
    sigmoidal func prevents highly enriched sgRNAs (e.g. jackpots) from having
    disproportionate influence vs. enriched but not jackpot sgRNAs
    input is the summed lfc (of sg1 and sg2), m, and theta
    theta controls the critical point (center), m controls steepness of function
    CLUMPS used m=3, theta=2, LSD1 used m=2, theta=3
    according to allison, m=3, theta=2 didn't work, so idk
    see the Gad Getz CLUMPS PNAS paper for reference
    """
    num = lfc**m
    denom = (lfc**m) + (theta**m)
    return num/denom

def calculate_pw_score(df_score, scores_col, tanh_a):
    """
    Calculate pairwise sums matrix for scores.
    """
    df_pws_scores = df_score.copy()

    # CALCULATE PAIRWISE SCORE OF EDITING DATA #
    pw_sum = df_pws_scores[scores_col].values[:, None] + df_pws_scores[scores_col].values[None, :] # PAIRWISE SUM MATRIX #
    df_pws_sum = pd.DataFrame(index=df_pws_scores['epegRNA_ID'],
                              columns=df_pws_scores['epegRNA_ID'], data=pw_sum)

    # CALC PARAMS FOR TANH OF Z-SCORE #
    upper_tri = np.where(np.triu(np.ones(df_pws_sum.shape), k=1).astype(bool),
                         df_pws_sum, np.nan) # GET UPPER TRIANGLE #
    pws_triu = pd.DataFrame(index=df_pws_sum.index,
                            columns=df_pws_sum.columns, data=upper_tri)

    flat_pws = pd.Series([y for x in pws_triu.columns for y in pws_triu[x]], name='sum_lfc').dropna() # FLATTEN FOR MEAN + STD CALC #
    df_pws = np.tanh(tanh_a * (df_pws_sum - flat_pws.mean()) / flat_pws.std()) # SCALED TANH FUNCTION #
    # COULD ALSO USE HILL FUNCTION HERE #
    return df_pws

### CALCULATE PWES

def calculate_pwes(df_gauss, df_pws, list_aas, pos_only):
    """
    Calculate PWES
    """
    df_pws.index, df_pws.columns = list_aas, list_aas
    # MAKE SURE THE INDEX OF df_pws AND df_gauss ARE THE SAME #
    df_pws = df_pws * df_gauss.loc[list_aas, list_aas].copy()

    if pos_only: df_pws[df_pws < 0] = 0 ###

    # SORT BY AA #
    df_pws_sort = df_pws.sort_index(axis=0)
    df_pws_sort = df_pws_sort.sort_index(axis=1)
    print("PWES calculated ...")

    return df_pws_sort, df_pws

def cluster_pwes(df_pws, df_score, list_aas, t, x_col):

    # Create dendrogram
    print("Starting linkage ...")
    df_clus = df_score.loc[df_score[x_col].isin(list_aas)].copy().reset_index()
    # df_clus['label_encoded'] = df_clus['label'].astype('category').cat.codes ###
    link = sp_cl.hierarchy.linkage(df_pws, method='ward', metric='euclidean', optimal_ordering=True)
    print("Linking complete ...")

    df_clus['cl_new'] = sp_cl.hierarchy.fcluster(link, t=t, criterion='distance')
    # Find number of clusters
    num_clus = sorted(df_clus['cl_new'].unique())[-1]
    print(f"Number of clusters: {num_clus}")

    return df_clus, link

def get_clus_aa(df_clus, x_col):
    """
    Given df_clus, generated from cluster_pwes, print the amino acids in each cluster.
    """
    try:
        df_clus['cl_new']
        df_clus[x_col]
    except KeyError:
        raise Exception('df_clus does not contain required columns (cl_new or aa_pos)')

    aas_dict = {}
    for clus in sorted(df_clus["cl_new"].unique()):
        df_subset = df_clus[df_clus["cl_new"] == clus][x_col]
        aas = sorted([x for x in list(set(df_subset))])
        print(f'Cluster {clus} amino acids: \n{aas}')
        aas_dict[clus] = aas
    return aas_dict

def shuffle_pwes(df_gauss, df_pws, list_aas, nrand=1000):
    """
    Randomize PWES Score
    """
    df_pws.index, df_pws.columns = list_aas, list_aas
    df_pws_list = []
    np.random.seed(0)

    for i in range(nrand):
        np.random.seed(i)
        permutation = np.random.permutation(df_pws.shape[0])
        df_shuffled_sample = df_pws.iloc[permutation, permutation] # SHUFFLE ROWS AND COLS #
        df_shuffled_sample.index, df_shuffled_sample.columns = list_aas, list_aas

        df_shuffled_sample *= df_gauss.loc[list_aas, list_aas].copy()
        df_pws_list.append(df_shuffled_sample)

    # COLLAPSE RANDOMIZATIONS INTO ONE MATRIX #
    result_df = sum(df_pws_list) / len(df_pws_list)
    # SORT BY AA #
    df_pws_sort = result_df.sort_index(axis=0).sort_index(axis=1)

    print("Randomized PWES calculated.")
    return df_pws_sort, result_df

### Iterative Clustering
def _safe_slug(x: Any, max_len: int = 120) -> str:
    """
    _safe_slug(): Make a filesystem-friendly folder name component.

    Parameters:
    x (Any): input value to convert to slug
    max_len (int, optional): maximum length of the resulting slug (Default: 120
    """
    s = str(x)
    s = s.strip()
    if s == "":
        s = "EMPTY"
    # Replace path separators and other weirdness
    s = s.replace(os.sep, "-").replace("/", "-")
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9._=-]+", "-", s)
    return s[:max_len]

def _mask_equals(series: pd.Series, value: Any) -> pd.Series:
    """
    _mask_equals(): Equality mask that treats NaN == NaN for grouping purposes.
    
    Parameters:
    series (pd.Series): pandas Series to compare
    value (Any): value to compare against (can be NaN)
    """
    if pd.isna(value):
        return series.isna()
    return series == value

### PWES Cluster Mapping
def pymol_script(
    dir: str,
    file: str,
    pdb_filename: str,
    residue_dict: dict,
    colors: list, 
    ):
    """
    pymol_script(): Generates a PyMOL script to display selected residues as separate objects.

    Parameters:
    dir (str): Path to save the generated PyMOL script.
    file (str): Name of the file to save the generated PyMOL script.
    pdb_filename (str): Path to the PDB file to load in PyMOL
    residue_dict (dict): Dictionary where keys are cluster identifiers and values are lists of residues (e.g., [('A', 123), ('B', 45)])
    colors (list): List of RGB tuples corresponding to each cluster for coloring the residues.
    """
    mkdir(dir)
    with open(os.path.join(dir, file), 'w') as pymol_script:
        pymol_script.write(f"load {pdb_filename}\n")
        pymol_script.write(f"show cartoon\n")
        pymol_script.write(f"color white\n")
        pymol_script.write(f"zoom\n")

        for i, (key, vals) in enumerate(residue_dict.items()):
            group_name = f"Cluster{key}"
            # res has a format of ('A', 123) #
            selection_str = " or ".join([f"resi {str(int(res[1:]))} and chain {res[0]} and name CA" for res in vals])

            color = colors[i]
            color_hex = tuple_to_hex(color)
            color_str = f'0x{color_hex}'
            print(color_str)
            
            pymol_script.write(f"select {selection_str}\n")
            pymol_script.write(f"create {group_name}, sele\n")
            pymol_script.write(f"show spheres, {group_name}\n")
            pymol_script.write(f"color {color_str}, {group_name}\n")
            pymol_script.write(f"disable {group_name}\n")

# Plot functions
def hist(df_clus: pd.DataFrame | str, cluster_col: str = "cl_new", line: float = None,
    file: str = None, dir: str = None, palette_or_cmap: str = 'turbo', alpha: float = 1.0, dodge: bool = False, jitter: bool = True, size: float = 5, edgecol: str = 'black', lw: int = 1, errorbar: str = 'sd', errwid: int = 1, errcap: float = 0.1,
    figsize: tuple = (5, 5), title: str = 'Number of Edits in Clusters', title_size: int = 12, title_weight: str = 'bold', title_font: str = 'Arial',
    x_axis: str = 'Cluster #', x_axis_size: int = 12, x_axis_weight: str = 'bold', x_axis_font: str = 'Arial', x_axis_scale: str = 'linear', x_axis_dims: tuple = (0, 0), x_axis_pad: int = None, x_ticks_size: int = 12, x_ticks_rot: int = 0, x_ticks_font: str = 'Arial', x_ticks: list = [],
    y_axis: str = 'Count', y_axis_size: int = 12, y_axis_weight: str = 'bold', y_axis_font: str = 'Arial', y_axis_scale: str = 'linear', y_axis_dims: tuple = (0, 0), y_axis_pad: int = None, y_ticks_size: int = 12, y_ticks_rot: int = 0, y_ticks_font: str = 'Arial', y_ticks: list = [],
    legend_title: str = '', legend_title_size: int = 12, legend_size: int = 12, legend_bbox_to_anchor: tuple = (1, 1), legend_loc: str = 'upper left', legend_items: tuple = (0, 0), legend_ncol: int = 1,
    legend_columnspacing: int=0, legend_handletextpad: float=0.5, legend_labelspacing: float=0.5, legend_borderpad: float=0.5, legend_handlelength: float=1, legend_size_html_multiplier: float=1.0,
    dpi: int = 0, transparent: bool = False, show: bool = True, space_capitalize: bool = True, **kwargs):
    """
    hist(): Generates a histogram of the number edits for each cluster.
    
    Parameters:
    df_clus (pd.DataFrame | str): DataFrame or path to a saved DataFrame containing clustering results, must have cluster_col column for cluster labels.
    cluster_col (str, optional): column name in df_clus with cluster labels (Default: "cl_new")
    line (float, optional): add horizontal line at y value or vertical line at x value
    file (str, optional): filename to save the plot (including the file extension)
    dir (str, optional): Directory to save the plot.
    palette_or_cmap (str, optional): Palette or colormap for the plot.
    alpha (float, optional): Alpha value for transparency of markers.
    edgecol (str, optional): Edge color for markers.
    figsize (tuple, optional): Figure size for the plot.
    title (str, optional): Title for the plot.
    title_size (int, optional): Size of the title font.
    title_weight (str, optional): Weight of the title font.
    title_font (str, optional): Font family of the title.
    x_axis (str, optional): Label for x-axis.
    x_axis_size (int, optional): Size of x-axis label font.
    x_axis_weight (str, optional): Weight of x-axis label font.
    x_axis_font (str, optional): Font family of x-axis label font.
    x_axis_scale (str, optional): Scale of x-axis (linear, log, etc.).
    x_axis_dims (tuple, optional): Dimensions of x-axis (min, max).
    x_axis_pad (int, optional): Padding for x-axis labels.
    x_ticks_size (int, optional): Size of x-tick labels font.
    x_ticks_rot (int, optional): Rotation angle of x-tick labels (in degrees).
    x_ticks_font (str, optional): Font family of x-tick labels font.
    x_ticks (list, optional): List of specific tick positions on the x-axis (if empty, default ticks are used).
    y_axis (str, optional): Label for y-axis.
    y_axis_size (int, optional): Size of y-axis label font.
    y_axis_weight (str, optional): Weight of y-axis label font.
    y_axis_font (str, optional): Font family of y-axis label font.
    y_axis_scale (str, optional): Scale of y-axis (linear, log, etc.).
    y_axis_dims (tuple, optional): Dimensions of y-axis (min, max).
    y_axis_pad (int, optional): Padding for y-axis labels.
    y_ticks_size (int, optional): Size of y-tick labels font.
    y_ticks_rot (int, optional): Rotation angle of y-tick labels (in degrees).
    y_ticks_font (str, optional): Font family of y-tick labels font.
    y_ticks (list, optional): List of specific tick positions on the y-axis (if empty, default ticks are used).
    legend_title (str, optional): Title for the legend.
    legend_title_size (int, optional): Size of the legend title font.
    legend_size (int, optional): Size of the legend labels font.
    legend_bbox_to_anchor (tuple, optional): Tuple specifying the position to anchor the legend (x, y).
    legend_loc (str, optional): Location of the legend (e.g., 'upper left', 'lower right').
    legend_items (tuple, optional): Tuple specifying the number of items in the legend (ncol, nrow).
    legend_ncol (int, optional): Number of columns in the legend.
    legend_columnspacing (float, optional): Spacing between columns in the legend.
    legend_handletextpad (float, optional): Padding between legend handles and text.
    legend_labelspacing (float, optional): Spacing between legend labels.
    legend_borderpad (float, optional): Padding between the legend border and its contents.
    legend_handlelength (float, optional): Length of the legend handles.
    legend_size_html_multiplier (float, optional): Multiplier for scaling legend font size when saving as HTML.
    dpi (int, optional): Dots per inch for saving the plot.
    transparent (bool, optional): Whether to save the plot with a transparent background.
    show (bool, optional): Whether to display the plot after creating it.
    space_capitalize (bool, optional): Whether to capitalize words in the title and axis labels.
    """
    if isinstance(df_clus, str):
        df_clus = io.get(df_clus)
    assert cluster_col in df_clus.columns, f"df_clus does not contain a '{cluster_col}' column"

    # Histogram of cluster counts
    clust_counts = df_clus[cluster_col].value_counts().sort_index().reset_index()
    clust_counts[cluster_col] = clust_counts[cluster_col].astype(str)  # Ensure cluster labels are strings for categorical plotting

    p.cat(graph='bar', df=clust_counts, x=cluster_col, y="count", line=line,
        file=f'{".".join(file.split(".")[:-1])}_hist.{file.split(".")[-1]}' if file else "hist.all", dir=dir,palette_or_cmap=palette_or_cmap,alpha=alpha,dodge=dodge,jitter=jitter,size=size,edgecol=edgecol,lw=lw,errorbar=errorbar,errwid=errwid,errcap=errcap,
        figsize=figsize,title=title,title_size=title_size,title_weight=title_weight,title_font=title_font,
        x_axis=x_axis,x_axis_size=x_axis_size,x_axis_weight=x_axis_weight,x_axis_font=x_axis_font,x_axis_scale=x_axis_scale,x_axis_dims=x_axis_dims,x_axis_pad=x_axis_pad,x_ticks_size=x_ticks_size,x_ticks_rot=x_ticks_rot,x_ticks_font=x_ticks_font,x_ticks=x_ticks,
        y_axis=y_axis,y_axis_size=y_axis_size,y_axis_weight=y_axis_weight,y_axis_font=y_axis_font,y_axis_scale=y_axis_scale,y_axis_dims=y_axis_dims,y_axis_pad=y_axis_pad,y_ticks_size=y_ticks_size,y_ticks_rot=y_ticks_rot,y_ticks_font=y_ticks_font,y_ticks=y_ticks,
        legend_title=legend_title,legend_title_size=legend_title_size,legend_size=legend_size,legend_bbox_to_anchor=legend_bbox_to_anchor,legend_loc=legend_loc,legend_items=legend_items,legend_ncol=legend_ncol,
        legend_columnspacing=legend_columnspacing,legend_handletextpad=legend_handletextpad,legend_labelspacing=legend_labelspacing,legend_borderpad=legend_borderpad,legend_handlelength=legend_handlelength,legend_size_html_multiplier=legend_size_html_multiplier,
        dpi=dpi,transparent=transparent,show=show,space_capitalize=space_capitalize,**kwargs)

def cat(df_clus: pd.DataFrame | str, cluster_col: str = "cl_new",scores_col: str = "log2(FC)", line: float = None,
    file: str = None, dir: str = None, palette_or_cmap: str = 'turbo', alpha: float = 1.0, dodge: bool = False, jitter: bool = True, size: float = 5, edgecol: str = 'black', lw: int = 1, errorbar: str = 'sd', errwid: int = 1, errcap: float = 0.1,
    figsize: tuple = (5, 5), title: str = '', title_size: int = 12, title_weight: str = 'bold', title_font: str = 'Arial',
    x_axis: str = '', x_axis_size: int = 12, x_axis_weight: str = 'bold', x_axis_font: str = 'Arial', x_axis_scale: str = 'linear', x_axis_dims: tuple = (0, 0), x_axis_pad: int = None, x_ticks_size: int = 12, x_ticks_rot: int = 0, x_ticks_font: str = 'Arial', x_ticks: list = [],
    y_axis: str = '', y_axis_size: int = 12, y_axis_weight: str = 'bold', y_axis_font: str = 'Arial', y_axis_scale: str = 'linear', y_axis_dims: tuple = (0, 0), y_axis_pad: int = None, y_ticks_size: int = 12, y_ticks_rot: int = 0, y_ticks_font: str = 'Arial', y_ticks: list = [],
    legend_title: str = '', legend_title_size: int = 12, legend_size: int = 12, legend_bbox_to_anchor: tuple = (1, 1), legend_loc: str = 'upper left', legend_items: tuple = (0, 0), legend_ncol: int = 1,
    legend_columnspacing: int=0, legend_handletextpad: float=0.5, legend_labelspacing: float=0.5, legend_borderpad: float=0.5, legend_handlelength: float=1, legend_size_html_multiplier: float=1.0,
    dpi: int = 0, transparent: bool = False, show: bool = True, space_capitalize: bool = True, **kwargs):
    """
    cat(): creates categorical graphs for PWES 3D clustering results
    
    Parameters:
    df_clus (pd.DataFrame | str): DataFrame or path to a saved DataFrame containing clustering results, must have cluster_col column for cluster labels and scores_col for values to plot.
    cluster_col (str, optional): column name in df_clus with cluster labels (Default: "cl_new")
    scores_col (str, optional): column name in df_clus with sgRNA scores (Default: "log2(FC)")
    line (float, optional): add horizontal line at y value or vertical line at x value
    file (str, optional): filename to save the plot (including the file extension)
    dir (str, optional): Directory to save the plot.
    palette_or_cmap (str, optional): Palette or colormap for the plot.
    alpha (float, optional): Alpha value for transparency of markers.
    edgecol (str, optional): Edge color for markers.
    figsize (tuple, optional): Figure size for the plot.
    title (str, optional): Title for the plot.
    title_size (int, optional): Size of the title font.
    title_weight (str, optional): Weight of the title font.
    title_font (str, optional): Font family of the title.
    x_axis (str, optional): Label for x-axis.
    x_axis_size (int, optional): Size of x-axis label font.
    x_axis_weight (str, optional): Weight of x-axis label font.
    x_axis_font (str, optional): Font family of x-axis label font.
    x_axis_scale (str, optional): Scale of x-axis (linear, log, etc.).
    x_axis_dims (tuple, optional): Dimensions of x-axis (min, max).
    x_axis_pad (int, optional): Padding for x-axis labels.
    x_ticks_size (int, optional): Size of x-tick labels font.
    x_ticks_rot (int, optional): Rotation angle of x-tick labels (in degrees).
    x_ticks_font (str, optional): Font family of x-tick labels font.
    x_ticks (list, optional): List of specific tick positions on the x-axis (if empty, default ticks are used).
    y_axis (str, optional): Label for y-axis.
    y_axis_size (int, optional): Size of y-axis label font.
    y_axis_weight (str, optional): Weight of y-axis label font.
    y_axis_font (str, optional): Font family of y-axis label font.
    y_axis_scale (str, optional): Scale of y-axis (linear, log, etc.).
    y_axis_dims (tuple, optional): Dimensions of y-axis (min, max).
    y_axis_pad (int, optional): Padding for y-axis labels.
    y_ticks_size (int, optional): Size of y-tick labels font.
    y_ticks_rot (int, optional): Rotation angle of y-tick labels (in degrees).
    y_ticks_font (str, optional): Font family of y-tick labels font.
    y_ticks (list, optional): List of specific tick positions on the y-axis (if empty, default ticks are used).
    legend_title (str, optional): Title for the legend.
    legend_title_size (int, optional): Size of the legend title font.
    legend_size (int, optional): Size of the legend labels font.
    legend_bbox_to_anchor (tuple, optional): Tuple specifying the position to anchor the legend (x, y).
    legend_loc (str, optional): Location of the legend (e.g., 'upper left', 'lower right').
    legend_items (tuple, optional): Tuple specifying the number of items in the legend (ncol, nrow).
    legend_ncol (int, optional): Number of columns in the legend.
    legend_columnspacing (float, optional): Spacing between columns in the legend.
    legend_handletextpad (float, optional): Padding between legend handles and text.
    legend_labelspacing (float, optional): Spacing between legend labels.
    legend_borderpad (float, optional): Padding between the legend border and its contents.
    legend_handlelength (float, optional): Length of the legend handles.
    legend_size_html_multiplier (float, optional): Multiplier for scaling legend font size when saving as HTML.
    dpi (int, optional): Dots per inch for saving the plot.
    transparent (bool, optional): Whether to save the plot with a transparent background.
    show (bool, optional): Whether to display the plot after creating it.
    space_capitalize (bool, optional): Whether to capitalize words in the title and axis labels.
    """
    if isinstance(df_clus, str):
        df_clus = io.get(df_clus)
    assert cluster_col in df_clus.columns, f"df_clus does not contain a '{cluster_col}' column"
    assert scores_col in df_clus.columns, f"{scores_col} not in df_clus"
    df_clus = df_clus.sort_values(cluster_col)  # Ensure clusters are sorted by their labels
    df_clus[cluster_col] = df_clus[cluster_col].astype(str)  # Ensure cluster labels are strings for categorical plotting
    x_axis = x_axis if x_axis else 'Cluster #'

    # Categorical bar chart of cluster counts
    p.cat(graph='swarm', df=df_clus, x=cluster_col, y=scores_col, line=line,
        file=f'{".".join(file.split(".")[:-1])}_swarm.{file.split(".")[-1]}' if file else f"swarm.all", dir=dir,palette_or_cmap=palette_or_cmap,alpha=alpha,dodge=dodge,jitter=jitter,size=size,edgecol=edgecol,lw=lw,errorbar=errorbar,errwid=errwid,errcap=errcap,
        figsize=figsize,title=title,title_size=title_size,title_weight=title_weight,title_font=title_font,
        x_axis=x_axis,x_axis_size=x_axis_size,x_axis_weight=x_axis_weight,x_axis_font=x_axis_font,x_axis_scale=x_axis_scale,x_axis_dims=x_axis_dims,x_axis_pad=x_axis_pad,x_ticks_size=x_ticks_size,x_ticks_rot=x_ticks_rot,x_ticks_font=x_ticks_font,x_ticks=x_ticks,
        y_axis=y_axis,y_axis_size=y_axis_size,y_axis_weight=y_axis_weight,y_axis_font=y_axis_font,y_axis_scale=y_axis_scale,y_axis_dims=y_axis_dims,y_axis_pad=y_axis_pad,y_ticks_size=y_ticks_size,y_ticks_rot=y_ticks_rot,y_ticks_font=y_ticks_font,y_ticks=y_ticks,
        legend_title=legend_title,legend_title_size=legend_title_size,legend_size=legend_size,legend_bbox_to_anchor=legend_bbox_to_anchor,legend_loc=legend_loc,legend_items=legend_items,legend_ncol=legend_ncol,
        legend_columnspacing=legend_columnspacing,legend_handletextpad=legend_handletextpad,legend_labelspacing=legend_labelspacing,legend_borderpad=legend_borderpad,legend_handlelength=legend_handlelength,legend_size_html_multiplier=legend_size_html_multiplier,
        dpi=dpi,transparent=transparent,show=show,space_capitalize=space_capitalize,**kwargs)
    
    p.cat(graph='strip', df=df_clus, x=cluster_col, y=scores_col, line=line,
        file=f'{".".join(file.split(".")[:-1])}_strip.{file.split(".")[-1]}' if file else f"strip.all", dir=dir,palette_or_cmap=palette_or_cmap,alpha=alpha,dodge=dodge,jitter=jitter,size=size,edgecol=edgecol,lw=lw,errorbar=errorbar,errwid=errwid,errcap=errcap,
        figsize=figsize,title=title,title_size=title_size,title_weight=title_weight,title_font=title_font,
        x_axis=x_axis,x_axis_size=x_axis_size,x_axis_weight=x_axis_weight,x_axis_font=x_axis_font,x_axis_scale=x_axis_scale,x_axis_dims=x_axis_dims,x_axis_pad=x_axis_pad,x_ticks_size=x_ticks_size,x_ticks_rot=x_ticks_rot,x_ticks_font=x_ticks_font,x_ticks=x_ticks,
        y_axis=y_axis,y_axis_size=y_axis_size,y_axis_weight=y_axis_weight,y_axis_font=y_axis_font,y_axis_scale=y_axis_scale,y_axis_dims=y_axis_dims,y_axis_pad=y_axis_pad,y_ticks_size=y_ticks_size,y_ticks_rot=y_ticks_rot,y_ticks_font=y_ticks_font,y_ticks=y_ticks,
        legend_title=legend_title,legend_title_size=legend_title_size,legend_size=legend_size,legend_bbox_to_anchor=legend_bbox_to_anchor,legend_loc=legend_loc,legend_items=legend_items,legend_ncol=legend_ncol,
        legend_columnspacing=legend_columnspacing,legend_handletextpad=legend_handletextpad,legend_labelspacing=legend_labelspacing,legend_borderpad=legend_borderpad,legend_handlelength=legend_handlelength,legend_size_html_multiplier=legend_size_html_multiplier,
        dpi=dpi,transparent=transparent,show=show,space_capitalize=space_capitalize,**kwargs)

    p.cat(graph='box', df=df_clus, x=cluster_col, y=scores_col, line=line,
        file=f'{".".join(file.split(".")[:-1])}_box.{file.split(".")[-1]}' if file else f"box.all", dir=dir,palette_or_cmap=palette_or_cmap,alpha=alpha,dodge=dodge,jitter=jitter,size=size,edgecol=edgecol,lw=lw,errorbar=errorbar,errwid=errwid,errcap=errcap,
        figsize=figsize,title=title,title_size=title_size,title_weight=title_weight,title_font=title_font,
        x_axis=x_axis,x_axis_size=x_axis_size,x_axis_weight=x_axis_weight,x_axis_font=x_axis_font,x_axis_scale=x_axis_scale,x_axis_dims=x_axis_dims,x_axis_pad=x_axis_pad,x_ticks_size=x_ticks_size,x_ticks_rot=x_ticks_rot,x_ticks_font=x_ticks_font,x_ticks=x_ticks,
        y_axis=y_axis,y_axis_size=y_axis_size,y_axis_weight=y_axis_weight,y_axis_font=y_axis_font,y_axis_scale=y_axis_scale,y_axis_dims=y_axis_dims,y_axis_pad=y_axis_pad,y_ticks_size=y_ticks_size,y_ticks_rot=y_ticks_rot,y_ticks_font=y_ticks_font,y_ticks=y_ticks,
        legend_title=legend_title,legend_title_size=legend_title_size,legend_size=legend_size,legend_bbox_to_anchor=legend_bbox_to_anchor,legend_loc=legend_loc,legend_items=legend_items,legend_ncol=legend_ncol,
        legend_columnspacing=legend_columnspacing,legend_handletextpad=legend_handletextpad,legend_labelspacing=legend_labelspacing,legend_borderpad=legend_borderpad,legend_handlelength=legend_handlelength,legend_size_html_multiplier=legend_size_html_multiplier,
        dpi=dpi,transparent=transparent,show=show,space_capitalize=space_capitalize,**kwargs)
    
    p.cat(graph='box_strip', df=df_clus, x=cluster_col, y=scores_col, line=line,
        file=f'{".".join(file.split(".")[:-1])}_box_strip.{file.split(".")[-1]}' if file else f"box_strip.all", dir=dir,palette_or_cmap=palette_or_cmap,alpha=alpha,dodge=dodge,jitter=jitter,size=size,edgecol=edgecol,lw=lw,errorbar=errorbar,errwid=errwid,errcap=errcap,
        figsize=figsize,title=title,title_size=title_size,title_weight=title_weight,title_font=title_font,
        x_axis=x_axis,x_axis_size=x_axis_size,x_axis_weight=x_axis_weight,x_axis_font=x_axis_font,x_axis_scale=x_axis_scale,x_axis_dims=x_axis_dims,x_axis_pad=x_axis_pad,x_ticks_size=x_ticks_size,x_ticks_rot=x_ticks_rot,x_ticks_font=x_ticks_font,x_ticks=x_ticks,
        y_axis=y_axis,y_axis_size=y_axis_size,y_axis_weight=y_axis_weight,y_axis_font=y_axis_font,y_axis_scale=y_axis_scale,y_axis_dims=y_axis_dims,y_axis_pad=y_axis_pad,y_ticks_size=y_ticks_size,y_ticks_rot=y_ticks_rot,y_ticks_font=y_ticks_font,y_ticks=y_ticks,
        legend_title=legend_title,legend_title_size=legend_title_size,legend_size=legend_size,legend_bbox_to_anchor=legend_bbox_to_anchor,legend_loc=legend_loc,legend_items=legend_items,legend_ncol=legend_ncol,
        legend_columnspacing=legend_columnspacing,legend_handletextpad=legend_handletextpad,legend_labelspacing=legend_labelspacing,legend_borderpad=legend_borderpad,legend_handlelength=legend_handlelength,legend_size_html_multiplier=legend_size_html_multiplier,
        dpi=dpi,transparent=transparent,show=show,space_capitalize=space_capitalize,**kwargs)

def torn(df: pd.DataFrame | str, df_clus: pd.DataFrame | str, cluster_col: str="cl_new", scores_col: str="log2(FC)", individual: bool=True, size: str | bool=None, size_dims: tuple=None, label: str='Edit', label_size: int=16,
        label_info: bool=True, aa_properties: bool | list=True, cBioPortal: str=None, only_clinical: bool=False, UniProt: str=None, PhosphoSitePlus: str=None, PDB_contacts: str=None, PDB_neighbors: str=None, ss_h: int=None, ss_y: int=None,
        file: str=None, dir: str=None, palette_or_cmap: str = 'turbo', edgecol: str='black', figsize=(5,5), title: str='', title_size: int = 12, title_weight: str='bold', title_font: str='Arial',
        x_axis: str='', x_axis_size: int=12, x_axis_weight: str='bold', x_axis_font: str='Arial', x_axis_dims: tuple=(0,0), x_axis_pad: int=None, x_ticks_size: int = 12, x_ticks_rot: int=0, x_ticks_font: str='Arial', x_ticks: list=[],
        y_axis: str='', y_axis_size: int=12, y_axis_weight: str='bold', y_axis_font: str='Arial', y_axis_dims: tuple=(0,0), y_axis_pad: int=None, y_ticks_size: int = 12, y_ticks_rot: int=0, y_ticks_font: str='Arial', y_ticks: list=[],
        legend_title: str='',legend_title_size: int=12, legend_size: int = 12, legend_bbox_to_anchor: tuple=(1,1), legend_loc: str='upper left', legend_ncol: int=1, 
        legend_columnspacing: int=-3, legend_handletextpad: float=0.5, legend_labelspacing: float=0.5, legend_borderpad: float=0.5, legend_handlelength: float=0.5, legend_size_html_multiplier: float=0.75,
        display_legend: bool=True, display_labels: bool=True, display_axis: bool=True, return_df: bool=True, dpi: int = 0, transparent: bool=False, show: bool=True, space_capitalize: bool=True,
        **kwargs) -> pd.DataFrame:
    ''' 
    torn(): creates tornado plots for visualizing PWES clusters
    
    Parameters:
    df (dataframe | str): pandas dataframe (or file path) from fastq.torn()
    df_clus (dataframe | str): pandas dataframe (or file path) from cluster_pwes()
    cluster_col (str, optional): column name in df & df_clus with cluster labels to plot (Default: 'cl_new')
    scores_col (str, optional): column name in df & df_clus with edit scores to plot (Default: 'log2(FC)')
    individual (bool, optional): show clusters on individual plot (otherwise, all clusters are shown on one plot)
    FC (str): fold change column name (y-axis)
    pval (str): p-value column name (size column if not specified)
    size (str | bool, optional): size column name (Default: pval; specify False for no size)
    size_dims (tuple, optional): (minimum,maximum) values in size column (Default: None)
    cluster_col (str, optional): cluster column name (Default: None)
    label (str, optional): label column name (Default: 'Edit'). Can't be None.
    label_size (int, optional): label font size (Default: 16)
    label_info (bool, optional): include additional info for labels if .html plot (Default: True)
    aa_properties (bool |list, optional): use aa_properties to format labels (Default: True). Options: True | False; ['hydrophobicity', 'polarity', 'charge', 'vdw_volume', 'pKa_C_term', 'pKa_N_term', 'pKa_side_chain']
    cBioPortal (str, optional): gene name (if saved to ~/.config/edms/cBioPortal_mutations) or file path for cBioPortal mutation data processed through edms.dat.cBioPortal.mutations()
    only_clinical (bool, optional): only show clinical mutations from cBioPortal (Default: False)
    UniProt (str, optional): UniProt accession (if saved to ~/.config/edms/UniProt) or file path for UniProt flat file. See edms.dat.uniprot.retrieve() or edms uniprot retrieve -h for more information.
    PhosphoSitePlus (str, optional): UniProt accession
    PDB_contacts (str, optional): PDB ID (if saved to ~/.config/edms/PDB) or file path for PDB structure file. See edms.dat.pdb.retrieve() or edms uniprot retrieve -h for more information.
    PDB_neighbors (str, optional): PDB ID (if saved to ~/.config/edms/PDB) or file path for PDB structure file. See edms.dat.pdb.retrieve() or edms uniprot retrieve -h for more information.
    ss_h (int, optional): height for secondary structure in the plot (Default: autogenerate)
    ss_y (int, optional): y position for secondary structure in the plot (Default: autogenerate)
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map (only used if individual parameter is specified)
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
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_axis_pad (int, optional): x-axis label padding
    x_ticks_size (int, optional): x-axis ticks font size
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-axis ticks font
    x_ticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_axis_pad (int, optional): y-axis label padding
    y_ticks_size (int, optional): y-axis ticks font size
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks_font (str, optional): y-axis ticks font
    y_ticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    legend_columnspacing (int, optional): space between columns (Default: -3; only for html plots)
    legend_handletextpad (float, optional): space between marker and text (Default: 0.5; only for html plots)
    legend_labelspacing (float, optional): vertical space between entries (Default: 0.5; only for html plots)
    legend_borderpad (float, optional): padding inside legend box (Default: 0.5; only for html plots)
    legend_handlelength (float, optional): marker length (Default: 0.5; only for html plots)
    legend_size_html_multiplier (float, optional): legend size multiplier for html plots (Default: 0.75)
    display_legend (bool, optional): display legend on plot (Default: True)
    display_labels (bool, optional): display labels for significant values (Default: True)
    display_axis (bool, optional): display x-axis line (Default: True)
    dpi (int, optional): figure dpi (Default: 1200 for non-HTML, 150 for HTML)
    transparent (bool, optional): save plot with transparent background (Default: False)
    return_df (bool, optional): return dataframe (Default: True)
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    
    Dependencies: os, matplotlib, seaborn, pandas, & edit_change()
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)
    if type(df_clus)==str:
        df_clus = io.get(pt=df_clus)
    
    # Organize data by conservation (changed to)
    stys_order = ['Conserved','Basic','Acidic','Polar','Nonpolar','Complex']
    mark_order = ['D','^','v','<','>','o']

    # Add label info: AA properties for conservation (change to); plus additional info from cBioPortal, UniProt, PhosphoSitePlus, PDB if specified
    df = add_label_info(df=df, label=label, label_size=label_size, label_info=label_info,
                        aa_properties=aa_properties, cBioPortal=cBioPortal, only_clinical=only_clinical, UniProt=UniProt, 
                        PhosphoSitePlus=PhosphoSitePlus, PDB_contacts=PDB_contacts, PDB_neighbors=PDB_neighbors)
    
    PDB_pt = None if PDB_contacts is None else PDB_contacts
    PDB_pt = PDB_pt if PDB_pt is not None else PDB_neighbors

    # Determine if we are saving to HTML (for interactive behavior)
    if file is not None:
        is_html = file.endswith('.html')

        if is_html == True:
            # Match title fontsize for html plots
            x_axis_size=title_size
            y_axis_size=title_size
            x_ticks_size=title_size
            y_ticks_size=title_size
            legend_title_size=title_size
            legend_size=title_size*legend_size_html_multiplier

            # Detailed labels for html plots
            if label_info == True:
                label = f'{label}_info'

    else:
        is_html = False

    # Organize data by 'pval' or specified 'size' column, typically input abundance
    sizes=(1,100)
    if size in [False,'False','false']: # No size
        size = None 

    else:
        if size is None: size = '-log10(pval)' # default to pval

        if size is not None and size in df.columns:
            # Filter by size dimensions
            if size_dims is not None: 
                df = df[(df[size]>=size_dims[0])&(df[size]<=size_dims[1])]

            # Shared size normalization across all scatter calls so marker areas are consistent
            size_norm = None
            _vmin, _vmax = None, None
            if display_legend:
                if size_dims is None:
                    _vmin = df[size].min()
                    _vmax = df[size].max()
                else:
                    _vmin = size_dims[0]
                    _vmax = size_dims[1]
                # Guard against degenerate case where all values are equal
                if _vmin == _vmax:
                    _vmax = _vmin + 1e-12
                size_norm = mcolors.Normalize(vmin=_vmin, vmax=_vmax)
    
    # Set dimensions
    if x_axis_dims==(0,0): x_axis_dims=(min(df[f'AA Number']),max(df[f'AA Number']))
    if y_axis_dims==(0,0): y_axis_dims=(min(df[scores_col]),max(df[scores_col]))

    if individual == True: # Plot all the clusters on individual plots
        for clus in df_clus[cluster_col].unique(): 
            clus_df = df_clus[df_clus[cluster_col]==clus]

            # Generate figure
            fig, ax = plt.subplots(figsize=figsize)

            # with data
            if display_legend==False: size=None
            stys='Change'
            sns.scatterplot(
                data=df,
                x='AA Number', y=scores_col,
                edgecolor=edgecol, color='gray', alpha=0.25,
                style=stys, style_order=stys_order if stys_order else None, markers=mark_order if mark_order else None,
                size=size if display_legend else None, sizes=sizes, size_norm=size_norm,
                legend=False, zorder = 1,
                ax=ax, **kwargs)
            sns.scatterplot(
                data=clus_df[clus_df[scores_col]<0],
                x='AA Number', y=scores_col,
                hue=scores_col, edgecolor=edgecol, palette='Blues_r', hue_norm=(df[scores_col].min(),0),
                style=stys, style_order=stys_order if stys_order else None, markers=mark_order if mark_order else None,
                size=size if display_legend else None, sizes=(50,50), size_norm=size_norm,
                legend=False, zorder = 2,
                ax=ax, **kwargs)
            sns.scatterplot(
                data=clus_df[clus_df[scores_col]>=0],
                x='AA Number', y=scores_col,
                hue=scores_col, edgecolor=edgecol, palette='Reds', hue_norm=(0,df[scores_col].max()),
                style=stys, style_order=stys_order if stys_order else None, markers=mark_order if mark_order else None,
                size=size if display_legend else None, sizes=(50,50), size_norm=size_norm, 
                legend=False, zorder = 3,
                ax=ax, **kwargs)
            
            # with x-axis line
            if display_axis == True:
                ax.plot([x_axis_dims[0], x_axis_dims[1]], [0,0], color='black', linestyle='-', linewidth=1)

            # with secondary structure
            if UniProt is not None:
                # Load UniProt flat file data
                try: # from filepath
                    UniProt_ss = uniprot.secondary_structure_from_flat_file(obj=UniProt)
                except:
                    try: # from config
                        for UniProt_file in os.listdir(os.path.expanduser('~/.config/edms/UniProt/')):
                            if UniProt.lower() in UniProt_file.lower():
                                UniProt_ss = uniprot.secondary_structure_from_flat_file(obj=f'{os.path.expanduser("~/.config/edms/UniProt")}/{UniProt_file}')
                                break
                    except:
                        raise FileNotFoundError(f"UniProt flat file not found: {UniProt}.\nPlease provide a valid filename or UniProt accession (if saved to {os.path.expanduser('~/.config/edms/UniProt/')}) or file path for UniProt flat file. See edms.dat.uniprot.retrieve_flat_file() or edms uniprot retrieve -h for more information.")
                
                # parameters for the secondary-structure "track"
                if ss_h is None:
                    ss_h  = 0.5 # height of secondary structure track
                if ss_y is None:
                    ss_y  = min(df[scores_col])-2*0.5 # position below min y value

                # helices
                helix_df = UniProt_ss[UniProt_ss["type"] == "α-helix"]
                if is_html:
                    helix_spans = [(row.start, row.end - row.start) for _, row in helix_df.iterrows()]
                    ax.broken_barh(helix_spans, (ss_y, ss_h), label="α-helix", facecolors='turquoise')
                else:
                    for xmin, xmax in t.zip_cols(df=helix_df, cols=['start','end']):
                        ax.axvspan(
                            xmin,
                            xmax,
                            facecolor='turquoise',
                        alpha=0.15,
                        edgecolor='none',
                        zorder=0)  # put behind scatter points

                # β-strands
                strand_df = UniProt_ss[UniProt_ss["type"] == "β-strand"]
                if is_html:
                    strand_spans = [(row.start, row.end - row.start) for _, row in strand_df.iterrows()]
                    ax.broken_barh(strand_spans, (ss_y, ss_h), label="β-strand", facecolors='gold')
                else:
                    for xmin, xmax in t.zip_cols(df=strand_df, cols=['start','end']):
                        ax.axvspan(
                            xmin,
                            xmax,
                            facecolor='gold',
                        alpha=0.15,
                        edgecolor='none',
                        zorder=0)  # put behind scatter points
                
            # with legend
            if display_legend == True:
                if size_norm is not None and size is not None: # Add consistent size legend with 5 representative values
                    if stys is not None and mark_order is not None and stys_order is not None: # Add stys legend too
                        legend_vals = np.linspace(_vmin, _vmax, len(mark_order))
                        for lv,mark,sty in zip(legend_vals,mark_order,stys_order):
                            plt.scatter([], [], s=np.interp(lv, [_vmin, _vmax], sizes), color='lightgray', label=f'{lv:.2g}; {sty}', marker=mark)
                        if is_html:
                            if legend_bbox_to_anchor == (1,1): legend_bbox_to_anchor = (-0.1,-0.2)
                            if legend_ncol == 1: legend_ncol = 3
                            plt.legend(title=legend_title if legend_title!='' else f'{size}; {stys}', 
                                        title_fontsize=legend_title_size, fontsize=legend_size,
                                        bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol,
                                        columnspacing=legend_columnspacing,    # space between columns
                                        handletextpad=legend_handletextpad,    # space between marker and text
                                        labelspacing=legend_labelspacing,      # vertical space between entries
                                        borderpad=legend_borderpad,            # padding inside legend box
                                        handlelength=legend_handlelength)      # marker length
                        else:
                            plt.legend(title=legend_title if legend_title!='' else f'{size}; {stys}', 
                                    title_fontsize=legend_title_size, fontsize=legend_size,
                                    bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol)
                    else:
                        legend_vals = np.linspace(_vmin, _vmax, 5)
                        for lv in legend_vals:
                            plt.scatter([], [], s=np.interp(lv, [_vmin, _vmax], sizes), color='lightgray', label=f'{lv:.2g}')
                        plt.legend(title=legend_title if legend_title!='' else size, 
                                    title_fontsize=legend_title_size, fontsize=legend_size,
                                    bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol)
            
            # with labels
            if display_labels == True:
                if is_html:
                    # For HTML, show labels interactively as tooltips instead of static text
                    pts = ax.scatter(
                        x=df_clus['AA Number'],
                        y=df_clus[scores_col],
                        s=20,
                        alpha=0
                    )
                    labels_list = df_clus[label].fillna("").astype(str).tolist()
                    tooltip = p.SafeHTMLTooltip(pts, labels_list)
                    clicker = p.ClickTooltip(pts, labels_list)
                    mpld3.plugins.connect(fig, tooltip, clicker)
                else:
                    # For static images, keep labels as always-visible text
                    for i, l in enumerate(df_clus[label]):
                        plt.text(
                            x=df_clus.iloc[i]['AA Number'],
                            y=df_clus.iloc[i][scores_col],
                            s=l
                        )
            
            # Set x axis
            if x_axis=='': x_axis='AA Number'
            plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight,fontfamily=x_axis_font, labelpad=x_axis_pad)
            if x_ticks==[]: 
                if x_ticks_rot==0: plt.xticks(rotation=x_ticks_rot,ha='center',va='top',fontfamily=x_ticks_font,fontsize=x_ticks_size)
                elif x_ticks_rot == 90: plt.xticks(rotation=x_ticks_rot,ha='right',va='center',fontfamily=x_ticks_font,fontsize=x_ticks_size)
                else: plt.xticks(rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font,fontsize=x_ticks_size)
            else: 
                if x_ticks_rot==0: plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot, ha='center',va='top',fontfamily=x_ticks_font,fontsize=x_ticks_size)
                elif x_ticks_rot == 90: plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot,ha='right',va='center',fontfamily=x_ticks_font,fontsize=x_ticks_size)
                else: plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font,fontsize=x_ticks_size)
            
            # Set y axis
            if y_axis=='': y_axis=f'{scores_col}'
            plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font, labelpad=y_axis_pad)

            if y_ticks==[]: plt.yticks(rotation=y_ticks_rot,fontfamily=y_ticks_font,fontsize=y_ticks_size)
            else: plt.yticks(ticks=y_ticks,labels=y_ticks,rotation=y_ticks_rot,fontfamily=y_ticks_font,fontsize=y_ticks_size)

            # Set title
            if title=='' and file is not None: 
                if space_capitalize: title=p.re_un_cap(".".join(file.split(".")[:-1]))
                else: title=".".join(file.split(".")[:-1])
            plt.title(title, fontsize=title_size, fontweight=title_weight, family=title_font)

            # Save & show fig; return dataframe
            p.save_fig(file=f"{'.'.join(file.split('.')[:-1])}_tornado_{clus}.{file.split('.')[-1]}" if file is not None else f"tornado.all", dir=dir, fig=fig, dpi=dpi, PDB_pt=PDB_pt, icon='tornado')
            if show:
                ext = file.split('.')[-1].lower() if file is not None else ''
                if ext not in ('html', 'json'):
                    plt.show()
                else:
                    mpld3.show(fig)
    
    else: # Only plot clusters above or below specified thresholds on the same plot
        clus_df = add_label_info(df=df_clus, label=label[:-5] if is_html else label, label_size=label_size, label_info=label_info,
                        aa_properties=aa_properties, cBioPortal=cBioPortal, only_clinical=only_clinical, UniProt=UniProt, 
                        PhosphoSitePlus=PhosphoSitePlus, PDB_contacts=PDB_contacts, PDB_neighbors=PDB_neighbors)

        # Generate figure
        fig, ax = plt.subplots(figsize=figsize)

        # with data
        if display_legend==False: size=None
        stys='Change'
        sns.scatterplot(
            data=df,
            x='AA Number', y=scores_col,
            edgecolor=edgecol, color='gray', alpha=0.25,
            style=stys, style_order=stys_order if stys_order else None, markers=mark_order if mark_order else None,
            size=size if display_legend else None, sizes=sizes, size_norm=size_norm,
            legend=False, zorder = 1,
            ax=ax, **kwargs)
        sns.scatterplot(
            data=clus_df,
            x='AA Number', y=scores_col,
            hue=cluster_col, edgecolor=edgecol, palette=palette_or_cmap, 
            style=stys, style_order=stys_order if stys_order else None, markers=mark_order if mark_order else None,
            size=size if display_legend else None, sizes=(50,50), size_norm=size_norm,
            legend=False, zorder = 2,
            ax=ax, **kwargs)
        
        # with x-axis line
        if display_axis == True:
            ax.plot([x_axis_dims[0], x_axis_dims[1]], [0,0], color='black', linestyle='-', linewidth=1)

        # with secondary structure
        if UniProt is not None:
            # Load UniProt flat file data
            try: # from filepath
                UniProt_ss = uniprot.secondary_structure_from_flat_file(obj=UniProt)
            except:
                try: # from config
                    for UniProt_file in os.listdir(os.path.expanduser('~/.config/edms/UniProt/')):
                        if UniProt.lower() in UniProt_file.lower():
                            UniProt_ss = uniprot.secondary_structure_from_flat_file(obj=f'{os.path.expanduser("~/.config/edms/UniProt")}/{UniProt_file}')
                            break
                except:
                    raise FileNotFoundError(f"UniProt flat file not found: {UniProt}.\nPlease provide a valid filename or UniProt accession (if saved to {os.path.expanduser('~/.config/edms/UniProt/')}) or file path for UniProt flat file. See edms.dat.uniprot.retrieve_flat_file() or edms uniprot retrieve -h for more information.")
            
            # parameters for the secondary-structure "track"
            if ss_h is None:
                ss_h  = 0.5 # height of secondary structure track
            if ss_y is None:
                ss_y  = min(df[scores_col])-2*0.5 # position below min y value

            # helices
            helix_df = UniProt_ss[UniProt_ss["type"] == "α-helix"]
            if is_html:
                helix_spans = [(row.start, row.end - row.start) for _, row in helix_df.iterrows()]
                ax.broken_barh(helix_spans, (ss_y, ss_h), label="α-helix", facecolors='turquoise')
            else:
                for xmin, xmax in t.zip_cols(df=helix_df, cols=['start','end']):
                    ax.axvspan(
                        xmin,
                        xmax,
                        facecolor='turquoise',
                    alpha=0.15,
                    edgecolor='none',
                    zorder=0)  # put behind scatter points

            # β-strands
            strand_df = UniProt_ss[UniProt_ss["type"] == "β-strand"]
            if is_html:
                strand_spans = [(row.start, row.end - row.start) for _, row in strand_df.iterrows()]
                ax.broken_barh(strand_spans, (ss_y, ss_h), label="β-strand", facecolors='gold')
            else:
                for xmin, xmax in t.zip_cols(df=strand_df, cols=['start','end']):
                    ax.axvspan(
                        xmin,
                        xmax,
                        facecolor='gold',
                    alpha=0.15,
                    edgecolor='none',
                    zorder=0)  # put behind scatter points
            
        # with legend
        if display_legend == True:
            if size_norm is not None and size is not None: # Add consistent size legend with 5 representative values
                if stys is not None and mark_order is not None and stys_order is not None: # Add stys legend too
                    legend_vals = np.linspace(_vmin, _vmax, len(mark_order))
                    for lv,mark,sty in zip(legend_vals,mark_order,stys_order): # Add legend entry for each size level with corresponding marker and stys
                        plt.scatter([], [], s=np.interp(lv, [_vmin, _vmax], sizes), color='lightgray', label=f'{lv:.2g}; {sty}', marker=mark)
                    #for clus in df_clus[cluster_col].unique(): # Add legend entry for each cluster
                        #plt.scatter([], [], color=palette_or_cmap[clus], label=f'Cluster {clus}', marker='o')
                    if is_html:
                        if legend_bbox_to_anchor == (1,1): legend_bbox_to_anchor = (-0.1,-0.2)
                        if legend_ncol == 1: legend_ncol = 3
                        plt.legend(title=legend_title if legend_title!='' else f'{size}; {stys}', 
                                    title_fontsize=legend_title_size, fontsize=legend_size,
                                    bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol,
                                    columnspacing=legend_columnspacing,    # space between columns
                                    handletextpad=legend_handletextpad,    # space between marker and text
                                    labelspacing=legend_labelspacing,      # vertical space between entries
                                    borderpad=legend_borderpad,            # padding inside legend box
                                    handlelength=legend_handlelength)      # marker length
                    else:
                        plt.legend(title=legend_title if legend_title!='' else f'{size}; {stys}', 
                                title_fontsize=legend_title_size, fontsize=legend_size,
                                bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol)
                else:
                    legend_vals = np.linspace(_vmin, _vmax, 5)
                    for lv in legend_vals:
                        plt.scatter([], [], s=np.interp(lv, [_vmin, _vmax], sizes), color='lightgray', label=f'{lv:.2g}')
                    plt.legend(title=legend_title if legend_title!='' else size, 
                                title_fontsize=legend_title_size, fontsize=legend_size,
                                bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol)
        
        # with labels
        if display_labels == True:
            if is_html:
                # For HTML, show labels interactively as tooltips instead of static text
                pts = ax.scatter(
                    x=clus_df['AA Number'],
                    y=clus_df[scores_col],
                    s=20,
                    alpha=0
                )
                labels_list = clus_df[label].fillna("").astype(str).tolist()
                tooltip = p.SafeHTMLTooltip(pts, labels_list)
                clicker = p.ClickTooltip(pts, labels_list)
                mpld3.plugins.connect(fig, tooltip, clicker)
            else:
                # For static images, keep labels as always-visible text
                for i, l in enumerate(clus_df[label]):
                    plt.text(
                        x=clus_df.iloc[i]['AA Number'],
                        y=clus_df.iloc[i][scores_col],
                        s=l
                    )
        
        # Set x axis
        if x_axis=='': x_axis='AA Number'
        plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight,fontfamily=x_axis_font, labelpad=x_axis_pad)
        if x_ticks==[]: 
            if x_ticks_rot==0: plt.xticks(rotation=x_ticks_rot,ha='center',va='top',fontfamily=x_ticks_font,fontsize=x_ticks_size)
            elif x_ticks_rot == 90: plt.xticks(rotation=x_ticks_rot,ha='right',va='center',fontfamily=x_ticks_font,fontsize=x_ticks_size)
            else: plt.xticks(rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font,fontsize=x_ticks_size)
        else: 
            if x_ticks_rot==0: plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot, ha='center',va='top',fontfamily=x_ticks_font,fontsize=x_ticks_size)
            elif x_ticks_rot == 90: plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot,ha='right',va='center',fontfamily=x_ticks_font,fontsize=x_ticks_size)
            else: plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font,fontsize=x_ticks_size)
        
        # Set y axis
        if y_axis=='': y_axis=f'{scores_col}'
        plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font, labelpad=y_axis_pad)

        if y_ticks==[]: plt.yticks(rotation=y_ticks_rot,fontfamily=y_ticks_font,fontsize=y_ticks_size)
        else: plt.yticks(ticks=y_ticks,labels=y_ticks,rotation=y_ticks_rot,fontfamily=y_ticks_font,fontsize=y_ticks_size)

        # Set title
        if title=='' and file is not None: 
            if space_capitalize: title=p.re_un_cap(".".join(file.split(".")[:-1]))
            else: title=".".join(file.split(".")[:-1])
        plt.title(title, fontsize=title_size, fontweight=title_weight, family=title_font)

        # Save & show fig; return dataframe
        p.save_fig(file=f"{'.'.join(file.split('.')[:-1])}_tornado.{file.split('.')[-1]}" if file is not None else f"tornado.all", dir=dir, fig=fig, dpi=dpi, transparent=transparent, PDB_pt=PDB_pt, icon='tornado')
        if show:
            ext = file.split('.')[-1].lower() if file is not None else ''
            if ext not in ('html', 'json'):
                plt.show()
            else:
                mpld3.show(fig)

    if return_df:
        return df

def heatmap(
    df_scaled: pd.DataFrame | str, mask_on: bool = True, blackout: list[str] = ['X'],
    file: str = None, dir: str = None, figsize: tuple = (4, 2.2),
    color: str='black', edgecol: str=None, lw: int=0, alpha: float=1,
    center: float=0, cmap: str="RdBu_r", sq: bool=True, vmin: float=-1, vmax: float=1, cbar: bool=False,
    spine_color: str = "black", spine_visible: bool = True, spine_wspace: float = 0.2, spine_hspace: float = 0.2,
    rasterize: bool=True, dpi: int=0, transparent: bool=False, show: bool=True,
    ):
    '''
    heatmap(): Generates a heatmap for the scaled PWES matrix with options for masking, blackout regions, and styling.

    Parameters:
    df_scaled (pd.DataFrame | str): Scaled DataFrame containing PWES values.
    mask_on (bool): Whether to mask the lower triangle of the heatmap.
    blackout (list[str]): List of prefixes for rows/columns to blackout (e.g., 'X' for non-protein residues).
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    figsize (tuple, optional): figure size per subplot
    color (str, optional): point color
    edgecol (str, optional): point edge color
    lw (int, optional): line width
    alpha (float, optional): transparency of blackout rectangles
    center (float, optional): value at which to center the colormap
    cmap (str, optional): colormap to use for heatmap
    sq (bool, optional): whether to make heatmap cells square
    vmin (float, optional): minimum value for colormap
    vmax (float, optional): maximum value for colormap
    cbar (bool, optional): whether to display colorbar
    spine_color (str, optional): color of plot spines
    spine_visible (bool, optional): whether to display plot spines
    spine_wspace (float, optional): width space for subplot grid (Default: 0.2)
    spine_hspace (float, optional): height space for subplot grid (Default: 0.2)
    rasterize (bool, optional): rasterize heatmap for smaller file size (Default: True)
    transparent (bool, optional): save plot with transparent background (Default: False)
    dpi (int, optional): figure dpi (Default: 1200 for non-HTML, 150 for HTML)
    show (bool, optional): show plot (Default: True)
    '''
    # Get dataframe if path is provided #
    if isinstance(df_scaled, str):
        df_scaled = io.get(df_scaled, index_col=0)

    # SETUP #
    mask = np.zeros_like(df_scaled)
    if mask_on: mask[np.tril_indices_from(mask, k=-1)] = True # MASK FOR LOWER TRIANGLE #

    # GENERATE HEATMAP #
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        data=df_scaled,
        ax=ax,
        mask=mask,
        center=center,
        square=sq,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        cbar=cbar,
        xticklabels=False,
        yticklabels=False,
        )

    for edge, spine in ax.spines.items():
        spine.set_visible(spine_visible)
        spine.set_color(spine_color)

    # SET TICKS AND LABELS #
    #ax.set_xlabel(x_axis, fontproperties={'family': x_axis_font, 'size': x_axis_size, 'weight': x_axis_weight})
    #ax.set_ylabel(y_axis, fontproperties={'family': y_axis_font, 'size': y_axis_size, 'weight': y_axis_weight})
    #ax.set_xticklabels(xticks)
    #ax.set_yticklabels(axis.yticks)
    # TICK LABELS #
    #for label in ax.get_xticklabels():
    #    label.set_fontproperties({'family': x_ticks_font, 'size': x_ticks_size, 'rotation': x_ticks_rot})
    #for label in ax.get_yticklabels(): 
    #    label.set_fontproperties({'family': y_ticks_font, 'size': y_ticks_size, 'rotation': y_ticks_rot})
    plt.subplots_adjust(wspace=spine_wspace, hspace=spine_hspace)

    # DRAW RECTANGLES #
    # HORIZONTAL #
    for i, idx in enumerate(df_scaled.index):
        for start in blackout:
            if idx.startswith(start):
                rect = patches.Rectangle(
                    (0, i), len(df_scaled.columns), 1,
                    linewidth=lw, edgecolor=edgecol, facecolor=color, alpha=alpha)
                ax.add_patch(rect)
    # VERTICAL #
    for j, col in enumerate(df_scaled.columns):
        for start in blackout:
            if col.startswith(start):
                rect = patches.Rectangle(
                    (j, 0), 1, len(df_scaled),
                    linewidth=lw, edgecolor=edgecol, facecolor=color, alpha=alpha)
                ax.add_patch(rect)

    # RASTERIZE #
    if rasterize:
        for coll in ax.collections:
            coll.set_rasterized(True)

    # SAVE & SHOW #
    plt.tight_layout()
    p.save_fig(file=file, dir=dir, fig=fig, dpi=dpi, transparent=transparent, icon='heat')
    if show:
        ext = file.split('.')[-1].lower() if file is not None else ''
        if ext not in ('html', 'json'):
            plt.show()
        else:
            mpld3.show(fig)
    plt.close(fig)

def heatmap_cbar(
    orientation: str='vertical', # otherwise horizontal
    file: str = None, dir: str = None, figsize: tuple = (.25, 2),
    cmap: str='RdBu_r', vmin: float=-1, vmax: float=1, lw: float = 0.5,
    tick_size: int=8, tick_width=0.5, label_font: str='Arial', label_size: int=6,
    dpi: int=0, transparent: bool=False, show: bool=True):
    '''
    heatmap_cbar(): Generates a standalone colorbar for the heatmap with specified orientation and styling.

    Parameters:
    orientation (str, optional): Orientation of the colorbar ('vertical' or 'horizontal').
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    figsize (tuple, optional): figure size per subplot
    cmap (str, optional): colormap to use for heatmap
    vmin (float, optional): minimum value for colormap
    vmax (float, optional): maximum value for colormap
    lw (float, optional): line width for colorbar outline (Default: 0.5)
    tick_size (int, optional): font size for colorbar ticks (Default: 8)
    tick_width (float, optional): tick width for colorbar (Default: 0.5)
    label_font (str, optional): label font (Default: Arial)
    label_size (int, optional): label font size (Default: 6)
    dpi (int, optional): figure dpi (Default: 1200 for non-HTML, 150 for HTML)
    transparent (bool, optional): save plot with transparent background (Default: False)
    show (bool, optional): show plot (Default: True)
    '''
    fig, ax = plt.subplots(figsize=figsize)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cbar = mpl.colorbar.ColorbarBase(
        ax, norm=norm, orientation=orientation,
        cmap=cmap)

    for label in (ax.get_yticklabels() if orientation == 'vertical' else ax.get_xticklabels()):
        label.set_fontproperties({'family': label_font, 'size': label_size})

    ax.tick_params(labelsize=tick_size, width=tick_width)
    cbar.outline.set_linewidth(lw)

    # SAVE & SHOW #
    p.save_fig(file=file, dir=dir, fig=fig, dpi=dpi, transparent=transparent, icon='heat')
    if show:
        ext = file.split('.')[-1].lower() if file is not None else ''
        if ext not in ('html', 'json'):
            plt.show()
        else:
            mpld3.show(fig)
    plt.close(fig)

def clustermap(
    df_scaled: pd.DataFrame | str, 
    link: np.ndarray | str, 
    df_clusters: pd.DataFrame | str,
    color_list: list=[],
    file: str = None, dir: str = None, figsize: tuple = (4, 2.2),
    center: float=0, cmap: str="RdBu_r", vmin: float=-1, vmax: float=1,
    line_color: str="k", line_style: str="--", lw: float=0.25,
    spine_color: str = "black", spine_visible: bool = True,
    rasterize: bool=True, dpi: int=0, transparent: bool=False, show: bool=True):
    '''
    pwes_clustermap(): Generates a clustered heatmap for the scaled PWES matrix with options for coloring clusters and styling.
    
    Parameters:
    df_scaled (pd.DataFrame | str): Scaled DataFrame containing PWES values.
    link (np.ndarray | str): Linkage matrix for hierarchical clustering.
    df_clusters (pd.DataFrame | str): DataFrame containing cluster assignments for rows.
    color_list (list): List of colors for clusters (optional).
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    figsize (tuple, optional): figure size per subplot
    center (float, optional): value at which to center the colormap
    cmap (str, optional): colormap to use for heatmap
    vmin (float, optional): minimum value for colormap
    vmax (float, optional): maximum value for colormap
    line_color (str, optional): line color
    line_style (str, optional): line style
    lw (float, optional): line width
    spine_color (str, optional): color of plot spines
    spine_visible (bool, optional): whether to display plot spines
    rasterize (bool, optional): rasterize heatmap for smaller file size (Default: True)
    dpi (int, optional): figure dpi (Default: 1200 for non-HTML, 150 for HTML)
    transparent (bool, optional): save plot with transparent background (Default: False)
    show (bool, optional): show plot (Default: True)
    '''
    # Get dataframes if paths are provided #
    if isinstance(df_scaled, str):
        df_scaled = io.get(df_scaled)
    if isinstance(df_clusters, str):
        df_clusters = io.get(df_clusters)
    
    # Get ndarray if linkage is provided as path #
    if isinstance(link, str):
        link = np.load(link)

    # SET COLORS FOR CLUSTERS #
    if color_list:
        df_colors = pd.DataFrame(index=df_scaled.index, columns=['Cluster'])
        df_colors['Colors'] = [color_list[i-1] for i in df_clusters['cl_new']]
        rcolor = df_colors['Colors'].copy()
    else: rcolor = None

    # CLUSTERMAP #
    fig = sns.clustermap(
        df_scaled,
        row_linkage=link, col_linkage=link,
        row_colors=rcolor,
        figsize=figsize,
        vmin=vmin,
        vmax=vmax,
        center=center, 
        cmap=cmap,
        xticklabels=False,
        yticklabels=False
        )
    ax = fig.ax_heatmap

    # AXIS SETTINGS #
    fig.ax_heatmap.collections[0].set_rasterized(rasterize)
    fig.cax.set_visible(False)
    ax.set_aspect('equal')
    fig.ax_col_dendrogram.set_visible(False)

    # DRAW LINES BETWEEN CLUSTERS #
    temp = 0
    clusterlist = df_clusters['cl_new'].value_counts().sort_index()
    for i in clusterlist.iloc[:-1]:
        temp += i
        ax.axhline(y=temp, color=line_color, linestyle=line_style, linewidth=lw)

    for edge, spine in ax.spines.items():
        spine.set_visible(spine_visible)
        spine.set_color(spine_color)

    # SAVE & SHOW #
    plt.tight_layout()
    p.save_fig(file=file, dir=dir, fig=fig, dpi=dpi, transparent=transparent, icon='heat')
    if show:
        ext = file.split('.')[-1].lower() if file is not None else ''
        if ext not in ('html', 'json'):
            plt.show()
        else:
            mpld3.show(fig)
    plt.close(fig)

# MAIN FUNCTION
def clustering(pdb_file: str, df: pd.DataFrame | str,
            x_col: str = "AA Number", scores_col: str = "log2(FC)",
            id_col: Optional[str] = None, iter_cols: Optional[Union[str, List[str]]] = None,
            gene_col: str = "", gene_map: Optional[Dict[str, str]] = None,
            tanh_a: float = 1, gauss_std: float = 16, dend_t: float = 10,
            aa_int: Optional[Tuple[int, int]] = None, plots: bool = True,
            out_prefix: Optional[str] = None, out_dir: Optional[str] = None,
            pos_only: bool = False, return_dfs: bool = False, pymol: bool = True, show: bool = False,
            hist_kwargs: Optional[Dict] = None, cat_kwargs: Optional[Dict] = None, torn_kwargs: Optional[Dict] = None,
            heatmap_kwargs: Optional[Dict] = None, heatmap_cbar_kwargs: Optional[Dict] = None, clustermap_kwargs: Optional[Dict] = None):
    """
    clustering(): Compute 3D PWES clustering analysis.

    Parameters:
    pdb_file (str): PDB id or path to PDB file for structural information. If PDB id is given, will attempt to retrieve from config or RCSB PDB.
    df (pd.DataFrame): pandas DataFrame with edit scores and AA positions
    x_col (str, optional): column name in df with AA positions (Default: "AA Number")
    scores_col (str, optional): column name in df with edit scores (Default: "log2(FC)")
    id_col (str, optional): column name in df with unique IDs for epegRNAs (if not using index)
    iter_cols (str | list[str], optional): If provided, run clustering separately for each value (or combination of values) in df[iter_cols]. Outputs are saved under: out_dir/<iter_cols>=<value>/  (or multi-col combo folder). The structure-derived distance matrices are computed once and reused.
    gene_col (str): column name in df with gene names (if multiple subunits)
    gene_map (dict, optional): dict mapping gene names to chain IDs in PDB file
    tanh_a (float, optional): parameter for tanh transformation of scores
    gauss_std (float, optional): standard deviation for Gaussian distance factor
    dend_t (float, optional): threshold for cutting dendrogram in clustering
    aa_int (tuple, optional): tuple of (aa_min, aa_max) defining the aas to calculate pdists. The default is None, which takes the min/max of df_centroids
    plots (bool, optional): if True, generate plots (Default: True)
    out_prefix (str, optional): prefix for output files (Default: timestamp)
    out_dir (str, optional): directory for output files (Default: "../out")
    pos_only (bool, optional): if True, only use positive scores for clustering (Default: False)
    return_dfs (bool, optional): if True, return all dataframes (Default: False)
    pymol (bool, optional): if True, generate PyMOL script(s) (Default: True)
    show (bool, optional): if True, display plots (Default: True)
    hist_kwargs (dict, optional): dict of keyword arguments to pass to hist() for histogram generation. If None, defaults will be used. (Default: None)
    cat_kwargs (dict, optional): dict of keyword arguments to pass to cat() for categorical plot generation. If None, defaults will be used. (Default: None)
    torn_kwargs (dict, optional): dict of keyword arguments to pass to torn() for tornado plot generation. If None, defaults will be used. (Default: None)
    heatmap_kwargs (dict, optional): dict of keyword arguments to pass to heatmap() for heatmap generation. If None, defaults will be used. (Default: None)
    heatmap_cbar_kwargs (dict, optional): dict of keyword arguments to pass to heatmap_cbar() for heatmap colorbar generation. If None, defaults will be used. (Default: None)
    clustmap_kwargs (dict, optional): dict of keyword arguments to pass to clustermap() for clustering generation. If None, defaults will be used. (Default: None)
    """
    gene_map = gene_map or {}

    # DEFAULTS + BASIC CHECKS #
    if out_dir is None:
        out_dir = "../out"
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    if out_prefix is None:
        out_prefix = datetime.now().strftime("%Y%m%d_%H%M%S")

    if isinstance(df, str): # Get df from path if needed
        df = io.get(df)
    assert x_col in df.columns, "Check [x_col] in df"
    assert scores_col in df.columns, "Check [scores_col] in df"

    # ---- CALCULATE SPATIAL COMPONENT FROM STRUCTURE (ONCE) ----
    if not Path(pdb_file).is_file():
        print(f"Could not find {pdb_file} at path. Attempting to retrieve from config...")
        try:
            for PDB_file in os.listdir(os.path.expanduser("~/.config/edms/PDB/")):
                if pdb_file.lower() in PDB_file.lower():
                    pdb_file = os.path.join(os.path.expanduser("~/.config/edms/PDB"), PDB_file)
                    break
        except Exception:
            print(f"Could not find {pdb_file} in config. Attempting to retrieve from RCSB PDB...")
            pdb.retrieve(id=pdb_file, suf=".pdb")
            pdb_file = os.path.join(os.path.expanduser("~/.config/edms/PDB"), pdb_file.lower() + ".pdb")

    df_centroids = process_pdb(pdb_file, gene_map.values())
    df_pwdist = get_pairwise_dist(df_centroids, gene_map.values(), aa_int)
    df_gauss = df_pwdist.apply(lambda x: gauss(x, gauss_std))
    del df_pwdist

    # -------------- inner runner: starts at ASSIGN IDs --------------
    def _run_one(df_sub: pd.DataFrame, out_dir_sub: str, out_prefix_sub: str):
        # ASSIGN IDs FOR EVERY INPUT #
        df_sub = df_sub.copy()

        if id_col is not None:
            assert id_col in df_sub.columns, "Check [id_col] in df"
            df_sub["epegRNA_ID"] = df_sub[id_col].astype(str)
        else:
            df_sub["epegRNA_ID"] = [
                "epegRNA" + num for num in map(str, list(range(df_sub.shape[0])))
            ]

        df_sub[x_col] = df_sub[x_col].astype(int)  # MUST BE INTEGER #

        # FOR ONLY ONE SUBUNIT #
        if (len(gene_map) == 0) or (not gene_col) or (gene_col == 0):
            print("No chain(s) indicated. Automatically assigning to structure ...")
            aa_in_structure = df_centroids.aa_num.tolist()
            df_sub = df_sub[df_sub[x_col].isin(aa_in_structure)]

            list_aas = df_sub[x_col]
            df_pws_score = calculate_pw_score(df_sub, scores_col, tanh_a)

        # FOR A COMPLEX #
        else:
            print("Chain(s) indicated. Mapping gene(s) to chains ...")
            assert gene_col in df_sub.columns, "Check [gene_col] is in df"

            aa_in_structure = df_centroids.label.tolist()
            df_sub = df_sub[df_sub[gene_col].isin(gene_map.keys())]  # FILTER OUT GENES NOT IN MAP #
            df_sub["label"] = (
                df_sub[gene_col].map(gene_map)
                + df_sub[x_col].astype(str).str.zfill(4)
            )
            df_sub = df_sub[df_sub["label"].isin(aa_in_structure)]

            list_aas = df_sub["label"]
            df_pws_score = calculate_pw_score(df_sub, scores_col, tanh_a)

        # CALCULATE PWES SCORE #
        df_pwes_sorted, df_pwes_unsorted = calculate_pwes(df_gauss, df_pws_score, list_aas, pos_only)

        # CLUSTER PWES #
        cluster_xcol = x_col if ((len(gene_map) == 0) or (not gene_col) or (gene_col == 0)) else "label"
        df_clus, link = cluster_pwes(
            df_pws=df_pwes_unsorted,
            df_score=df_sub,
            list_aas=list_aas,
            t=dend_t,
            x_col=cluster_xcol,
        )

        # SAVE CLUSTERS + DATAFRAMES #
        Path(out_dir_sub).mkdir(parents=True, exist_ok=True)

        aas_dict = get_clus_aa(df_clus, cluster_xcol)
        with open(os.path.join(out_dir_sub, f"{out_prefix_sub}_aas_dict.pickle"), "wb") as f:
            pickle.dump(aas_dict, f)

        out = os.path.join(out_dir_sub, out_prefix_sub)
        df_gauss.to_csv(f"{out}_df_gauss.csv")                # shared, but saved per group for convenience
        df_pws_score.to_csv(f"{out}_df_pairwise.csv")
        df_pwes_sorted.to_csv(f"{out}_df_pwes_sorted.csv")
        df_pwes_unsorted.to_csv(f"{out}_df_pwes_unsorted.csv")
        df_clus.to_csv(f"{out}_df_clus.csv")
        np.savetxt(f"{out}_linkage.npy", np.asarray(link, dtype=float))

        # Generate plots if requested
        if plots:
            if hist_kwargs is None:
                hist(df_clus=df_clus, dir=out_dir_sub, title=f'{out_prefix.replace("_", " ")}\n{out_dir_sub.split("/")[-1].replace("__", "; ").replace("_", " ")}', file = f"{out_prefix_sub}.all", show=show)
            else:
                hist(df_clus=df_clus, dir=out_dir_sub, title=f'{out_prefix.replace("_", " ")}\n{out_dir_sub.split("/")[-1].replace("__", "; ").replace("_", " ")}', file = f"{out_prefix_sub}.all", show=show, **hist_kwargs)
            if cat_kwargs is None:
                cat(df_clus=df_clus, scores_col=scores_col, dir=out_dir_sub, title=f'{out_prefix.replace("_", " ")}\n{out_dir_sub.split("/")[-1].replace("__", "; ").replace("_", " ")}', file = f"{out_prefix_sub}.all", show=show)
            else:
                cat(df_clus=df_clus, scores_col=scores_col, dir=out_dir_sub, title=f'{out_prefix.replace("_", " ")}\n{out_dir_sub.split("/")[-1].replace("__", "; ").replace("_", " ")}', file = f"{out_prefix_sub}.all", show=show, **cat_kwargs)
            if torn_kwargs is None:
                torn(df=df, df_clus=df_clus, dir=out_dir_sub, title=f'{out_prefix.replace("_", " ")}\n{out_dir_sub.split("/")[-1].replace("__", "; ").replace("_", " ")}', file = f"{out_prefix_sub}.all", show=show)
            else:
                torn(df=df, df_clus=df_clus, dir=out_dir_sub, title=f'{out_prefix.replace("_", " ")}\n{out_dir_sub.split("/")[-1].replace("__", "; ").replace("_", " ")}', file = f"{out_prefix_sub}.all", show=show, **torn_kwargs)
            if heatmap_kwargs is None:
                heatmap(df_scaled=df_pwes_sorted, dir=out_dir_sub, file = f"{out_prefix_sub}_heatmap.all", show=show)
            else:
                heatmap(df_scaled=df_pwes_sorted, dir=out_dir_sub, file = f"{out_prefix_sub}_heatmap.all", show=show, **heatmap_kwargs)
            if heatmap_cbar_kwargs is None:
                heatmap_cbar(dir=out_dir_sub, file = f"{out_prefix_sub}_heatmap_colorbar.all", show=show)
            else:
                heatmap_cbar(dir=out_dir_sub, file = f"{out_prefix_sub}_heatmap_colorbar.all", show=show, **heatmap_cbar_kwargs)
            if clustermap_kwargs is None:
                clustermap(df_scaled=df_pwes_sorted, link=link, df_clusters=df_clus, dir=out_dir_sub, file = f"{out_prefix_sub}_clustermap.all", show=show)
            else:
                clustermap(df_scaled=df_pwes_sorted, link=link, df_clusters=df_clus, dir=out_dir_sub, file = f"{out_prefix_sub}_clustermap.all", show=show, **clustermap_kwargs)
        
        ### FIGURE THIS OUT ###
        #if pymol:
        #    pymol_script(dir=out_dir_sub, file=f"{out_prefix_sub}_pymol.txt",
        #                 pdb_filename=pdb_file, residue_dict=aas_dict, colors="tab20")

        if return_dfs:
            return df_gauss, df_pwes_sorted, df_pwes_unsorted, df_clus, aas_dict, link
        return None

    # ---------------- iterative mode ----------------
    if iter_cols is not None:
        cols = [iter_cols] if isinstance(iter_cols, str) else list(iter_cols)
        for c in cols:
            assert c in df.columns, f"iter_cols column not found: {c}"

        # value_counts on single col or combination (MultiIndex)
        vc = df[cols].value_counts(dropna=False)

        results = {}  # group_key -> outputs (only if return_dfs=True)

        for key_vals, _count in vc.items():
            print(f"\nRunning group: {dict(zip(cols, key_vals))} ({_count} rows)")

            # Normalize key_vals to tuple aligned with cols
            if len(cols) == 1:
                key_tuple = (key_vals,)
            else:
                key_tuple = tuple(key_vals)

            # Build mask
            mask = pd.Series(True, index=df.index)
            for c, v in zip(cols, key_tuple):
                mask &= _mask_equals(df[c], v)

            df_sub = df.loc[mask].copy()
            if df_sub.empty:
                continue

            # Subfolder name
            if len(cols) == 1:
                folder = f"{cols[0]}={_safe_slug(key_tuple[0])}"
            else:
                parts = [f"{c}={_safe_slug(v)}" for c, v in zip(cols, key_tuple)]
                folder = "__".join(parts)

            out_dir_sub = os.path.join(out_dir, folder)
            try:
                res = _run_one(df_sub, out_dir_sub=out_dir_sub, out_prefix_sub=out_prefix)
            except Exception as e:
                print(f"Error running _run_one for group {folder}: {e}")
                res = None
            if return_dfs:
                results[key_tuple if len(cols) > 1 else key_tuple[0]] = res

        if return_dfs:
            return results
        return None

    # ---------------- normal (non-iterative) mode ----------------
    return _run_one(df, out_dir_sub=out_dir, out_prefix_sub=out_prefix)