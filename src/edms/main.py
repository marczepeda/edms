#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
'''
Module: main.py
Author: Marc Zepeda
Created: 2025-04-12
Description: Endogenous Deep Mutational Scans

Usage:
[Show post-install notice only once per version]
- try: importlib.metadata to get current version
- _user_config_dir(appname): Get cross-platform user config directory
- _maybe_show_first_run_notice(appname): Show a post-install notice only once per version

[Supporting argument methods]
- parse_tuple_int(arg): Parse a string argument into a tuple of integers
- parse_tuple_float(arg): Parse a string argument into a tuple of floats

[Plot subparser methods]
- add_common_plot_scat_args(subparser): Add common arguments for scatter plot related graphs
- add_common_plot_cat_args(subparser): Add common arguments for category dependent graphs
- add_common_plot_dist_args(subparser): Add common arguments for distribution graphs
- add_common_plot_heat_args(subparser): Add common arguments for heatmap graphs
- add_common_plot_stack_args(subparser): Add common arguments for stacked bar plot
- add_common_plot_vol_args(subparser): Add common arguments for volcano plot

[Main method]
- main(): Endogenous Deep Mutational Scans
'''
# Import packages
from __future__ import annotations # NEEDS TO BE FIRST LINE
import json, os, sys
from pathlib import Path
import argparse
import argcomplete
import ast
import datetime

from . import config
from . import utils

from .gen import plot as p
from .gen import stat as st
from .gen import io
from .gen import cli
from .gen import tidy as t

from .dat import cosmic as co 
from .dat import cvar

from .bio import ngs
from .bio import transfect as tf
from .bio import qPCR
from .bio import clone as cl
from .bio import fastq as fq
from .bio import pe
from .bio import sanger

# Show post-install notice only once per version
try:
    # Py3.8+: importlib.metadata in stdlib
    from importlib.metadata import version, PackageNotFoundError
except Exception:  # pragma: no cover
    version = lambda _: "0.0.0"  # fallback

def _user_config_dir(appname: str) -> Path:
    """
    _user_config_dir(): Cross-platform config dir without extra deps.
    
    Parameters:
    appname (str): Name of the application (used to create app-specific subdirectory)
    """
    if os.name == "nt":  # Windows
        base = Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming"))
    else:  # POSIX
        base = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
    return base / appname

def _maybe_show_first_run_notice(appname: str = "edms") -> None:
    """
    _maybe_show_first_run_notice(): Show a post-install notice only once per version.
    
    Parameters:
    appname (str): Name of the application (used to create app-specific subdirectory)
    """
    # Allow users/CI to suppress the message
    if os.environ.get("EDMS_NO_FIRST_RUN") == "1":
        return

    cfg_dir = _user_config_dir(appname)
    state_file = cfg_dir / "state.json"
    cfg_dir.mkdir(parents=True, exist_ok=True)

    # Determine current installed version (best effort)
    try:
        cur_ver = version("edms")
    except Exception:
        cur_ver = None

    state = {}
    if state_file.exists():
        try:
            state = json.loads(state_file.read_text() or "{}")
        except json.JSONDecodeError:
            state = {}

    last_shown_for = state.get("first_run_notice_shown_for")

    # Show only if never shown for this version (or ever)
    if cur_ver and last_shown_for == cur_ver:
        return

    # ---- Your message goes here ----
    print(t.color("\nrecommended post-install:\n\n  1) edms autocomplete\n     Follow CLI instructions or open new terminal\n\n  2) conda create -n umi_tools umi_tools cutadapt samtools bowtie2 fgbio\n     Required for edms fastq {extract_umis,trim_motifs,make_sams,make_bams,bam_umi_tags,group_umis,consensus_umis,bam_to_fastq}\n", "yellow"),
        file=sys.stderr,
    )
    # --------------------------------

    state["first_run_notice_shown_for"] = cur_ver or "unknown"
    try:
        state_file.write_text(json.dumps(state, indent=2))
    except OSError:
        # Non-fatal: if we can't write, just skip persisting
        pass

# Supporting argument methods
def parse_tuple_int(arg):
    '''
    parse_tuple_int(arg): Parse a string argument into a tuple of integers
    '''
    try:
        return tuple(map(int, arg.split(',')))
    except:
        raise argparse.ArgumentTypeError(f"{arg} must be formatted as 'int,int,...'")
    
def parse_tuple_float(arg):
    '''
    parse_tuple_float(arg): Parse a string argument into a tuple of floats
    '''
    try:
        return tuple(map(float, arg.split(',')))
    except:
        raise argparse.ArgumentTypeError(f"{arg} must be formatted as 'float,float,...'")

# Plot subparser methods
def add_common_plot_scat_args(subparser):
    '''
    add_common_plot_scat_args(subparser): Add common arguments for scatter plot related graphs
    '''
    # scat(): Required arguments
    subparser.add_argument("--df", help="Input dataframe file path", type=str, required=True)
    subparser.add_argument("--x", help="X-axis column", type=str, required=True)
    subparser.add_argument("--y", help="Y-axis column", type=str, required=True)

    # Optional core arguments
    subparser.add_argument("--cols", type=str, help="Color column name")
    subparser.add_argument("--cols_ord", nargs="+", help="Column order (list of values)")
    subparser.add_argument("--cols_exclude", nargs="+", help="Columns to exclude from coloring")
    subparser.add_argument("--stys", type=str, help="Style column name")

    subparser.add_argument("--dir", help="Output directory path", type=str, default='./out')
    subparser.add_argument("--file", help="Output file name", type=str, required=False, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_scat.png')
    subparser.add_argument("--palette_or_cmap", type=str, default="colorblind", help="Seaborn palette or matplotlib colormap")
    subparser.add_argument("--edgecol", type=str, default="black", help="Edge color for scatter points")

    # Figure appearance
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(10,6), help="Figure size as a tuple: width,height")
    subparser.add_argument("--title", type=str, default="", help="Plot title")
    subparser.add_argument("--title_size", type=int, default=18, help="Plot title font size")
    subparser.add_argument("--title_weight", type=str, default="bold", help="Plot title font weight (e.g., bold, normal)")
    subparser.add_argument("--title_font", type=str, default="Arial", help="Font family for the title")

    # X-axis settings
    subparser.add_argument("--x_axis", type=str, default="", help="X-axis label")
    subparser.add_argument("--x_axis_size", type=int, default=12, help="X-axis label font size")
    subparser.add_argument("--x_axis_weight", type=str, default="bold", help="X-axis label font weight (e.g., bold)")
    subparser.add_argument("--x_axis_font", type=str, default="Arial", help="X-axis label font family")
    subparser.add_argument("--x_axis_scale", type=str, default="linear", help="X-axis scale: linear, log, etc.")
    subparser.add_argument("--x_axis_dims", type=parse_tuple_int, default=(0,0), help="X-axis range as a tuple: start,end")
    subparser.add_argument("--x_ticks_rot", type=int, default=0, help="Rotation angle of X-axis tick labels")
    subparser.add_argument("--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")
    subparser.add_argument("--x_ticks", nargs="+", help="Specific tick values for X-axis")

    # Y-axis settings
    subparser.add_argument("--y_axis", type=str, default="", help="Y-axis label")
    subparser.add_argument("--y_axis_size", type=int, default=12, help="Y-axis label font size")
    subparser.add_argument("--y_axis_weight", type=str, default="bold", help="Y-axis label font weight")
    subparser.add_argument("--y_axis_font", type=str, default="Arial", help="Y-axis label font family")
    subparser.add_argument("--y_axis_scale", type=str, default="linear", help="Y-axis scale: linear, log, etc.")
    subparser.add_argument("--y_axis_dims", type=parse_tuple_int, default=(0,0), help="Y-axis range as a tuple: start,end")
    subparser.add_argument("--y_ticks_rot", type=int, default=0, help="Rotation angle of Y-axis tick labels")
    subparser.add_argument("--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")
    subparser.add_argument("--y_ticks", nargs="+", help="Specific tick values for Y-axis")

    # Legend settings
    subparser.add_argument("--legend_title", type=str, default="", help="Legend title")
    subparser.add_argument("--legend_title_size", type=int, default=12, help="Legend title font size")
    subparser.add_argument("--legend_size", type=int, default=9, help="Legend font size")
    subparser.add_argument("--legend_bbox_to_anchor", type=parse_tuple_float, default=(1,1), help="Bounding box anchor position for legend")
    subparser.add_argument("--legend_loc", type=str, default="upper left", help="Location of the legend in the plot")
    subparser.add_argument("--legend_items", type=parse_tuple_int, default=(0,0), help="Legend item count as a tuple (used for layout)")
    subparser.add_argument("--legend_ncol", type=int, default=1, help="Number of columns in legend")

    # Display and formatting
    subparser.add_argument("--show", action="store_true", help="Show the plot", default=False)
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize label/legend strings and replace underscores with spaces")

def add_common_plot_cat_args(subparser):
    '''
    add_common_plot_cat_args(subparser): Add common arguments for category dependent graphs
    '''
    # cat(): Required arguments
    subparser.add_argument("--df", help="Input dataframe file path", type=str, required=True)

    # Optional core arguments
    subparser.add_argument("--x", help="X-axis column name", type=str, default="")
    subparser.add_argument("--y", help="Y-axis column name", type=str, default="")
    subparser.add_argument("--cols", type=str, help="Color column name for grouping")
    subparser.add_argument("--cols_ord", nargs="+", help="Color column values order")
    subparser.add_argument("--cols_exclude", nargs="+", help="Color column values to exclude")

    subparser.add_argument("--file", type=str, help="Output filename", default='./out')
    subparser.add_argument("--dir", type=str, help="Output directory", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_cat.png')
    subparser.add_argument("--palette_or_cmap", type=str, default="colorblind", help="Seaborn color palette or matplotlib colormap")
    subparser.add_argument("--edgecol", type=str, default="black", help="Edge color for markers")

    # Error bar and style options
    subparser.add_argument("--lw", type=int, default=1, help="Line width for plot edges")
    subparser.add_argument("--errorbar", type=str, default="sd", help="Error bar type: sd (standard deviation), etc.")
    subparser.add_argument("--errwid", type=float, default=1, help="Width of the error bars")
    subparser.add_argument("--errcap", type=float, default=0.1, help="Cap size on error bars")

    # Figure appearance
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(10,6), help="Figure size formatted 'width,height'")
    subparser.add_argument("--title", type=str, default="", help="Plot title text")
    subparser.add_argument("--title_size", type=int, default=18, help="Font size of the plot title")
    subparser.add_argument("--title_weight", type=str, default="bold", help="Font weight of the plot title (e.g., bold, normal)")
    subparser.add_argument("--title_font", type=str, default="Arial", help="Font family for the plot title")

    # X-axis settings
    subparser.add_argument("--x_axis", type=str, default="", help="X-axis label text")
    subparser.add_argument("--x_axis_size", type=int, default=12, help="Font size for the X-axis label")
    subparser.add_argument("--x_axis_weight", type=str, default="bold", help="Font weight for the X-axis label")
    subparser.add_argument("--x_axis_font", type=str, default="Arial", help="Font family for the X-axis label")
    subparser.add_argument("--x_axis_scale", type=str, default="linear", help="Scale of X-axis (e.g., linear, log)")
    subparser.add_argument("--x_axis_dims", type=parse_tuple_float, default=(0, 0), help="X-axis range as tuple: start,end")
    subparser.add_argument("--x_ticks_rot", type=int, default=0, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")
    subparser.add_argument("--x_ticks", nargs="+", help="Explicit tick values for X-axis")

    # Y-axis settings
    subparser.add_argument("--y_axis", type=str, default="", help="Y-axis label text")
    subparser.add_argument("--y_axis_size", type=int, default=12, help="Font size for the Y-axis label")
    subparser.add_argument("--y_axis_weight", type=str, default="bold", help="Font weight for the Y-axis label")
    subparser.add_argument("--y_axis_font", type=str, default="Arial", help="Font family for the Y-axis label")
    subparser.add_argument("--y_axis_scale", type=str, default="linear", help="Scale of Y-axis (e.g., linear, log)")
    subparser.add_argument("--y_axis_dims", type=parse_tuple_float, default=(0, 0), help="Y-axis range as tuple: start,end")
    subparser.add_argument("--y_ticks_rot", type=int, default=0, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")
    subparser.add_argument("--y_ticks", nargs="+", help="Explicit tick values for Y-axis")

    # Legend settings
    subparser.add_argument("--legend_title", type=str, default="", help="Title for the legend")
    subparser.add_argument("--legend_title_size", type=int, default=12, help="Font size for the legend title")
    subparser.add_argument("--legend_size", type=int, default=9, help="Font size for legend items")
    subparser.add_argument("--legend_bbox_to_anchor", type=parse_tuple_float, default=(1, 1), help="Anchor position of the legend bounding box")
    subparser.add_argument("--legend_loc", type=str, default="upper left", help="Location of the legend on the plot")
    subparser.add_argument("--legend_items", type=parse_tuple_int, default=(0, 0), help="Tuple for legend item layout")
    subparser.add_argument("--legend_ncol", type=int, default=1, help="Number of columns in the legend")

    # Display and formatting
    subparser.add_argument("--show", action="store_true", help="Show the plot in a window", default=False)
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize labels and replace underscores with spaces")

def add_common_plot_dist_args(subparser):
    '''
    add_common_plot_dist_args(subparser): Add common arguments for distribution graphs
    '''
    # dist(): Required argument
    subparser.add_argument("--df", help="Input dataframe file path", type=str, required=True)
    subparser.add_argument("--x", type=str, help="X-axis column name", required=True)

    # File output
    subparser.add_argument("--dir", type=str, help="Output directory", default='./out')
    subparser.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_dist.png')

    # Optional core arguments
    subparser.add_argument("--cols", type=str, help="Color column name for grouping")
    subparser.add_argument("--cols_ord", nargs="+", help="Custom order for color column values")
    subparser.add_argument("--cols_exclude", nargs="+", help="Color column values to exclude")

    # Plot customization
    subparser.add_argument("--bins", type=int, default=40, help="Number of bins for histogram")
    subparser.add_argument("--log10_low", type=int, default=0, help="Log10 scale lower bound (e.g., 1 = 10^1 = 10)")
    subparser.add_argument("--palette_or_cmap", type=str, default="colorblind", help="Seaborn color palette or matplotlib colormap")
    subparser.add_argument("--edgecol", type=str, default="black", help="Edge color of histogram bars")
    subparser.add_argument("--lw", type=int, default=1, help="Line width for edges")
    subparser.add_argument("--ht", type=float, default=1.5, help="Height of the plot")
    subparser.add_argument("--asp", type=int, default=5, help="Aspect ratio of the plot")
    subparser.add_argument("--tp", type=float, default=0.8, help="Top padding space")
    subparser.add_argument("--hs", type=int, default=0, help="Horizontal spacing between plots (if faceted)")
    subparser.add_argument("--des", action="store_true", help="Remove plot spines (despine)", default=False)

    # Figure appearance
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(10,6), help="Figure size formatted as 'width,height'")
    subparser.add_argument("--title", type=str, default="", help="Plot title text")
    subparser.add_argument("--title_size", type=int, default=18, help="Plot title font size")
    subparser.add_argument("--title_weight", type=str, default="bold", help="Plot title font weight (e.g., bold, normal)")
    subparser.add_argument("--title_font", type=str, default="Arial", help="Font family for the plot title")

    # X-axis
    subparser.add_argument("--x_axis", type=str, default="", help="Label for the X-axis")
    subparser.add_argument("--x_axis_size", type=int, default=12, help="Font size for X-axis label")
    subparser.add_argument("--x_axis_weight", type=str, default="bold", help="Font weight for X-axis label")
    subparser.add_argument("--x_axis_font", type=str, default="Arial", help="Font family for X-axis label")
    subparser.add_argument("--x_axis_scale", type=str, default="linear", help="X-axis scale (e.g., linear, log)")
    subparser.add_argument("--x_axis_dims", type=parse_tuple_float, default=(0, 0), help="X-axis range as tuple: start,end")
    subparser.add_argument("--x_ticks_rot", type=int, default=0, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")
    subparser.add_argument("--x_ticks", nargs="+", help="Explicit tick values for X-axis")

    # Y-axis
    subparser.add_argument("--y_axis", type=str, default="", help="Label for the Y-axis")
    subparser.add_argument("--y_axis_size", type=int, default=12, help="Font size for Y-axis label")
    subparser.add_argument("--y_axis_weight", type=str, default="bold", help="Font weight for Y-axis label")
    subparser.add_argument("--y_axis_font", type=str, default="Arial", help="Font family for Y-axis label")
    subparser.add_argument("--y_axis_scale", type=str, default="linear", help="Y-axis scale (e.g., linear, log)")
    subparser.add_argument("--y_axis_dims", type=parse_tuple_float, default=(0, 0), help="Y-axis range as tuple: start,end")
    subparser.add_argument("--y_ticks_rot", type=int, default=0, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")
    subparser.add_argument("--y_ticks", nargs="+", help="Explicit tick values for Y-axis")

    # Legend
    subparser.add_argument("--legend_title", type=str, default="", help="Title text for the legend")
    subparser.add_argument("--legend_title_size", type=int, default=12, help="Font size of the legend title")
    subparser.add_argument("--legend_size", type=int, default=9, help="Font size for legend items")
    subparser.add_argument("--legend_bbox_to_anchor", type=parse_tuple_float, default=(1, 1), help="Legend bbox anchor position")
    subparser.add_argument("--legend_loc", type=str, default="upper left", help="Legend location on the plot")
    subparser.add_argument("--legend_items", type=parse_tuple_int, default=(0, 0), help="Tuple for legend layout items")
    subparser.add_argument("--legend_ncol", type=int, default=1, help="Number of columns in the legend")

    # Final display
    subparser.add_argument("--show", action="store_true", help="Show the plot in an interactive window", default=False)
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize and space legend/label values", default=False)

def add_common_plot_heat_args(subparser):
    '''
    add_common_plot_heat_args(subparser): Add common arguments for heatmap graphs
    '''
    # Required arguments
    subparser.add_argument("--df", help="Input dataframe file path", type=str, required=True)

    # Optional arguments
    subparser.add_argument("--x", type=str, help="X-axis column name to pivot tidy-formatted dataframe into matrix format")
    subparser.add_argument("--y", type=str, help="Y-axis column name to pivot tidy-formatted dataframe into matrix format")
    subparser.add_argument("--vars", type=str, help="Variable column name to split tidy-formatted dataframe into a dictionary of pivoted dataframes")
    subparser.add_argument("--vals", type=str, help="Value column name to populate pivoted dataframes")
    subparser.add_argument("--vals_dims", type=parse_tuple_float, help="Value column limits formatted as 'vmin,vmax'")

    subparser.add_argument("--dir", type=str, help="Output directory path", default='./out')
    subparser.add_argument("--file", type=str, help="Output filename", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_heat.png')
    subparser.add_argument("--edgecol", type=str, default="black", help="Color of cell edges")
    subparser.add_argument("--lw", type=int, default=1, help="Line width for cell borders")

    subparser.add_argument("--annot", action="store_true", help="Display cell values as annotations", default=False)
    subparser.add_argument("--cmap", type=str, default="Reds", help="Matplotlib colormap to use for heatmap")
    subparser.add_argument("--sq", action="store_true", help="Use square aspect ratio for cells", default=False)
    subparser.add_argument("--cbar", action="store_true", help="Display colorbar", default=False)

    # Title and size
    subparser.add_argument("--title", type=str, default="", help="Plot title")
    subparser.add_argument("--title_size", type=int, default=18, help="Font size of the title")
    subparser.add_argument("--title_weight", type=str, default="bold", help="Font weight of the title (e.g., bold, normal)")
    subparser.add_argument("--title_font", type=str, default="Arial", help="Font family for the title")
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(10,6), help="Figure size formatted as 'width,height'")

    # X-axis
    subparser.add_argument("--x_axis", type=str, default="", help="X-axis label")
    subparser.add_argument("--x_axis_size", type=int, default=12, help="Font size for X-axis label")
    subparser.add_argument("--x_axis_weight", type=str, default="bold", help="Font weight for X-axis label")
    subparser.add_argument("--x_axis_font", type=str, default="Arial", help="Font family for X-axis label")
    subparser.add_argument("--x_ticks_rot", type=int, default=0, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")

    # Y-axis
    subparser.add_argument("--y_axis", type=str, default="", help="Y-axis label")
    subparser.add_argument("--y_axis_size", type=int, default=12, help="Font size for Y-axis label")
    subparser.add_argument("--y_axis_weight", type=str, default="bold", help="Font weight for Y-axis label")
    subparser.add_argument("--y_axis_font", type=str, default="Arial", help="Font family for Y-axis label")
    subparser.add_argument("--y_ticks_rot", type=int, default=0, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")

    # Final display
    subparser.add_argument("--show", action="store_true", help="Show the plot in an interactive window", default=False)
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize and space labels/legend values", default=False)

def add_common_plot_stack_args(subparser):
    '''
    add_common_plot_stack_args(subparser): Add common arguments for stacked bar plot
    '''
    # Required arguments
    subparser.add_argument("--df", type=str, help="Input dataframe file path", required=True)
    subparser.add_argument("--x", type=str, help="X-axis column name")
    subparser.add_argument("--y", type=str, help="Y-axis column name")
    subparser.add_argument("--cols", type=str, help="Color column name for stacking")

    # Optional parameters
    subparser.add_argument("--dir", type=str, help="Output directory path", default='./out')
    subparser.add_argument("--file", type=str, help="Output filename", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_stack.png')
    
    subparser.add_argument("--cutoff", type=float, default=0, help="Minimum value cutoff for inclusion in stack")
    subparser.add_argument("--cols_ord", nargs="+", help="Order of values in the color column")
    subparser.add_argument("--x_ord", nargs="+", help="Custom order of X-axis categories")
    subparser.add_argument("--palette_or_cmap", type=str, default="Set2", help="Seaborn palette or Matplotlib colormap for stacked bars")
    subparser.add_argument("--errcap", type=int, default=4, help="Width of error bar caps")
    subparser.add_argument("--vertical", action="store_true", help="Stack bars vertically (default True)", default=False)

    # Figure & layout
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(10,6), help="Figure size formatted as 'width,height'")
    subparser.add_argument("--title", type=str, default="", help="Plot title")
    subparser.add_argument("--title_size", type=int, default=18, help="Font size of the title")
    subparser.add_argument("--title_weight", type=str, default="bold", help="Font weight of the title (e.g., bold, normal)")
    subparser.add_argument("--title_font", type=str, default="Arial", help="Font family for the title")

    # X-axis formatting
    subparser.add_argument("--x_axis", type=str, default="", help="X-axis label")
    subparser.add_argument("--x_axis_size", type=int, default=12, help="Font size for X-axis label")
    subparser.add_argument("--x_axis_weight", type=str, default="bold", help="Font weight for X-axis label")
    subparser.add_argument("--x_axis_font", type=str, default="Arial", help="Font family for X-axis label")
    subparser.add_argument("--x_ticks_rot", type=int, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")

    # Y-axis formatting
    subparser.add_argument("--y_axis", type=str, default="", help="Y-axis label")
    subparser.add_argument("--y_axis_size", type=int, default=12, help="Font size for Y-axis label")
    subparser.add_argument("--y_axis_weight", type=str, default="bold", help="Font weight for Y-axis label")
    subparser.add_argument("--y_axis_font", type=str, default="Arial", help="Font family for Y-axis label")
    subparser.add_argument("--y_axis_dims", type=parse_tuple_float, default=(0,0), help="Y-axis range as tuple: start,end")
    subparser.add_argument("--y_ticks_rot", type=int, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")

    # Legend options
    subparser.add_argument("--legend_title", type=str, default="", help="Legend title text")
    subparser.add_argument("--legend_title_size", type=int, default=12, help="Font size of the legend title")
    subparser.add_argument("--legend_size", type=int, default=12, help="Font size for legend items")
    subparser.add_argument("--legend_bbox_to_anchor", type=parse_tuple_float, default=(1, 1), help="Anchor position for the legend bounding box")
    subparser.add_argument("--legend_loc", type=str, default="upper left", help="Legend location on the plot")
    subparser.add_argument("--legend_ncol", type=int, default=1, help="Number of columns in the legend")

    # Display and formatting
    subparser.add_argument("--show", action="store_true", help="Show the plot in an interactive window", default=False)
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize and space legend/label values", default=False)

def add_common_plot_vol_args(subparser):
    '''
    add_common_plot_vol_args(subparser): Add common arguments for volcano plot
    '''
    # Required arguments
    subparser.add_argument("--df", type=str, help="Input dataframe file path")
    subparser.add_argument("--x", type=str, help="X-axis column name (e.g., fold change)")
    subparser.add_argument("--y", type=str, help="Y-axis column name (e.g., p-value)")

    # Optional data columns
    subparser.add_argument("--stys", type=str, help="Style column name for custom markers")
    subparser.add_argument("--size", type=str, help="Column name used to scale point sizes")
    subparser.add_argument("--size_dims", type=parse_tuple_float, help="Size range for points formatted as min,max")
    subparser.add_argument("--label", type=str, help="Column containing text labels for points")

    # Thresholds
    subparser.add_argument("--FC_threshold", type=float, default=2, help="Fold change threshold for significance")
    subparser.add_argument("--pval_threshold", type=float, default=0.05, help="P-value threshold for significance")

    # Output
    subparser.add_argument("--dir", type=str, help="Output directory path", default='./out')
    subparser.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_vol.png')

    # Aesthetics
    subparser.add_argument("--color", type=str, default="lightgray", help="Color for non-significant points")
    subparser.add_argument("--alpha", type=float, default=0.5, help="Transparency for non-significant points")
    subparser.add_argument("--edgecol", type=str, default="black", help="Edge color of points")
    subparser.add_argument("--vertical", action="store_true", help="Use vertical layout for plot", default=False)

    # Figure setup
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(10,6), help="Figure size formatted as 'width,height'")
    subparser.add_argument("--title", type=str, default="", help="Plot title")
    subparser.add_argument("--title_size", type=int, default=18, help="Font size for plot title")
    subparser.add_argument("--title_weight", type=str, default="bold", help="Font weight for plot title (e.g., bold, normal)")
    subparser.add_argument("--title_font", type=str, default="Arial", help="Font family for plot title")

    # X-axis settings
    subparser.add_argument("--x_axis", type=str, default="", help="Label for the X-axis")
    subparser.add_argument("--x_axis_size", type=int, default=12, help="Font size for X-axis label")
    subparser.add_argument("--x_axis_weight", type=str, default="bold", help="Font weight for X-axis label")
    subparser.add_argument("--x_axis_font", type=str, default="Arial", help="Font family for X-axis label")
    subparser.add_argument("--x_axis_dims", type=parse_tuple_float, default=(0, 0), help="X-axis range formatted as min,max")
    subparser.add_argument("--x_ticks_rot", type=int, default=0, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")
    subparser.add_argument("--x_ticks", nargs="+", help="Custom tick values for X-axis")

    # Y-axis settings
    subparser.add_argument("--y_axis", type=str, default="", help="Label for the Y-axis")
    subparser.add_argument("--y_axis_size", type=int, default=12, help="Font size for Y-axis label")
    subparser.add_argument("--y_axis_weight", type=str, default="bold", help="Font weight for Y-axis label")
    subparser.add_argument("--y_axis_font", type=str, default="Arial", help="Font family for Y-axis label")
    subparser.add_argument("--y_axis_dims", type=parse_tuple_float, default=(0, 0), help="Y-axis range formatted as min,max")
    subparser.add_argument("--y_ticks_rot", type=int, default=0, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")
    subparser.add_argument("--y_ticks", nargs="+", help="Custom tick values for Y-axis")

    # Legend
    subparser.add_argument("--legend_title", type=str, default="", help="Title for the legend")
    subparser.add_argument("--legend_title_size", type=int, default=12, help="Font size for legend title")
    subparser.add_argument("--legend_size", type=int, default=9, help="Font size for legend items")
    subparser.add_argument("--legend_bbox_to_anchor", type=parse_tuple_float, default=(1, 1), help="Bounding box anchor for legend")
    subparser.add_argument("--legend_loc", type=str, default="upper left", help="Legend location on the plot")
    subparser.add_argument("--legend_items", type=parse_tuple_int, default=(0, 0), help="Tuple defining legend layout")
    subparser.add_argument("--legend_ncol", type=int, default=1, help="Number of columns in the legend")

    # Boolean switches
    subparser.add_argument("--display_size", action="store_true", help="Show point sizes as annotations", default=False)
    subparser.add_argument("--display_labels", action="store_true", help="Display text labels for selected points", default=False)
    subparser.add_argument("--return_df", action="store_true", help="Return annotated DataFrame after plotting", default=False)
    subparser.add_argument("--show", action="store_true", help="Show the plot in an interactive window", default=False)
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize and space labels/legend items", default=False)

# Main method
def main():
    '''
    main(): Endogenous Deep Mutational Scans
    '''
    print(t.color("project: Endogenous Deep Mutational Scans (EDMS)", "green"))
    _maybe_show_first_run_notice()

    # Add parser and subparsers
    parser = argparse.ArgumentParser(prog="edms", description="Endogenous Deep Mutational Scans (EDMS)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(dest="command") # dest="command" required for autocomplete

    '''
    edms.utils:
    - run_bundled_script(): Run a bundled script from the package script resources.
    '''
    subparsers.add_parser("autocomplete", help="Enable edms autocomplete") # Run autocomplete script (args.command == "autocomplete")

    '''
    edms.config:
    - get_info: Retrieve information based on id
    - set_info: Set information based on id
    '''
    parser_config = subparsers.add_parser("config", help="Configuration")
    subparsers_config = parser_config.add_subparsers()

    # get_info(): Retrieve information based on id
    # set_info(): Set information based on id
    parser_config_get_info = subparsers_config.add_parser("get", help="Retrieve information based on id")
    parser_config_set_info = subparsers_config.add_parser("set", help="Set information based on id")
    
    # get_info() arguments
    parser_config_get_info.add_argument("--id", type=str, help="Identifier from/for configuration file")
    
    # set_info() arguments
    parser_config_set_info.add_argument("--id", type=str, help="Identifier from/for configuration file", required=True)
    parser_config_set_info.add_argument("--info", type=ast.literal_eval, help="Information for configuration file (str or dict)", required=True)
    
    # default functions
    parser_config_get_info.set_defaults(func=config.get_info)
    parser_config_set_info.set_defaults(func=config.set_info)

    '''
    edms.gen.plot:
    - scat(): creates scatter plot related graphs
    - cat(): creates category dependent graphs
    - dist(): creates distribution graphs
    - heat(): creates heatmap graphs
    - stack(): creates stacked bar plot
    - vol(): creates volcano plot
    '''
    parser_plot = subparsers.add_parser("plot", help="Generate scatter, category, distribution, heatmap, stacked bar, and volcano plots")
    subparsers_plot = parser_plot.add_subparsers(dest="typ")

    # scat(): Creates scatter plot related graphs (scat, line, line_scat)
    parser_plot_type_scat = subparsers_plot.add_parser("scat", help="Create scatter plot")
    parser_plot_type_line = subparsers_plot.add_parser("line", help="Create line plot")
    parser_plot_type_line_scat = subparsers_plot.add_parser("line_scat", help="Create scatter + line plot")

    for parser_plot_scat in [parser_plot_type_scat, parser_plot_type_line, parser_plot_type_line_scat]:
        add_common_plot_scat_args(parser_plot_scat)
        parser_plot_scat.set_defaults(func=p.scat)

    # cat(): Creates category dependent graphs (bar, box, violin, swarm, strip, point, count, bar_swarm, box_swarm, violin_swarm)
    parser_plot_type_bar = subparsers_plot.add_parser("bar", help="Create bar plot")
    parser_plot_type_box = subparsers_plot.add_parser("box", help="Create box plot")
    parser_plot_type_violin = subparsers_plot.add_parser("violin", help="Create violin plot")
    parser_plot_type_swarm = subparsers_plot.add_parser("swarm", help="Create swarm plot")
    parser_plot_type_strip = subparsers_plot.add_parser("strip", help="Create strip plot")
    parser_plot_type_point = subparsers_plot.add_parser("point", help="Create point plot")
    parser_plot_type_count = subparsers_plot.add_parser("count", help="Create count plot")
    parser_plot_type_bar_swarm = subparsers_plot.add_parser("bar_swarm", help="Create bar + swarm plot")
    parser_plot_type_box_swarm = subparsers_plot.add_parser("box_swarm", help="Create box + swarm plot")
    parser_plot_type_violin_swarm = subparsers_plot.add_parser("violin_swarm", help="Create violin + swarm plot")

    for parser_plot_cat in [parser_plot_type_bar, parser_plot_type_box, parser_plot_type_violin, parser_plot_type_swarm, parser_plot_type_strip, parser_plot_type_point, parser_plot_type_count, parser_plot_type_bar_swarm, parser_plot_type_box_swarm, parser_plot_type_violin_swarm]:
        add_common_plot_cat_args(parser_plot_cat)
        parser_plot_cat.set_defaults(func=p.cat)

    # dist(): Creates distribution graphs (hist, kde, hist_kde, rid)
    parser_plot_type_hist = subparsers_plot.add_parser("hist", help="Create histogram plot")
    parser_plot_type_kde = subparsers_plot.add_parser("kde", help="Create density plot")
    parser_plot_type_hist_kde = subparsers_plot.add_parser("hist_kde", help="Create histogram + density plot")
    parser_plot_type_rid = subparsers_plot.add_parser("rid", help="Create ridge plot")

    for parser_plot_dist in [parser_plot_type_hist, parser_plot_type_kde, parser_plot_type_hist_kde, parser_plot_type_rid]:
        add_common_plot_dist_args(parser_plot_dist)
        parser_plot_cat.set_defaults(func=p.dist)

    # heat(): Creates heatmap graphs
    parser_plot_type_heat = subparsers_plot.add_parser("heat", help="Create heatmap plot")
    add_common_plot_heat_args(parser_plot_type_heat)
    parser_plot_type_heat.set_defaults(func=p.heat)
    
    # stack(): Creates stacked bar plot
    parser_plot_type_stack = subparsers_plot.add_parser("stack", help="Create stacked bar plot")
    add_common_plot_stack_args(parser_plot_type_stack)
    parser_plot_type_stack.set_defaults(func=p.stack)

    # vol(): Creates volcano plot
    parser_plot_type_vol = subparsers_plot.add_parser("vol", help="Create volcano plot")
    add_common_plot_vol_args(parser_plot_type_vol)
    parser_plot_type_vol.set_defaults(func=p.vol)

    '''
    edms.gen.stat:
    - describe(): returns descriptive statistics for numerical columns in a DataFrame
    - difference(): computes the appropriate statistical test(s) and returns the p-value(s)
    - correlation(): returns a correlation matrix
    - compare(): computes FC, pval, and log transformations relative to a specified condition
    '''
    parser_stat = subparsers.add_parser("stat", help="Statistics")
    subparsers_stat = parser_stat.add_subparsers()
    
    # describe(): returns descriptive statistics for numerical columns in a DataFrame
    parser_stat_describe = subparsers_stat.add_parser("describe", help="Compute descriptive statistics")

    parser_stat_describe.add_argument("--df", type=str, help="Input file path", required=True)

    parser_stat_describe.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_stat_describe.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_descriptive.csv')
    
    parser_stat_describe.add_argument("--cols", nargs="+", help="List of numerical columns to describe")
    parser_stat_describe.add_argument("--group", type=str, help="Column name to group by")
    
    parser_stat_describe.set_defaults(func=st.describe)

    # difference(): computes the appropriate statistical test(s) and returns the p-value(s)
    parser_stat_difference = subparsers_stat.add_parser("difference", help="Compute statistical difference between groups")

    parser_stat_difference.add_argument("--df", type=str, help="Input file path",required=True)
    parser_stat_difference.add_argument("--data_col", type=str, help="Name of column containing numerical data",required=True)
    parser_stat_difference.add_argument("--compare_col", type=str, help="Name of column used for grouping/comparisons",required=True)
    parser_stat_difference.add_argument("--compare", nargs="+", help="List of groups to compare (e.g. A B)",required=True)

    parser_stat_difference.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_stat_difference.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_difference.csv')

    parser_stat_difference.add_argument("--same", action="store_true", help="Same subjects (paired test)")
    parser_stat_difference.add_argument("--para", action="store_true", help="Use parametric test (Default: True)")
    parser_stat_difference.add_argument("--alpha", type=float, default=0.05, help="Significance level (Default: 0.05)")
    parser_stat_difference.add_argument("--within_cols", nargs="+", help="Columns for repeated measures (used if same=True and para=True)")
    parser_stat_difference.add_argument("--method", type=str, default="holm", help="Correction method for multiple comparisons")

    parser_stat_difference.set_defaults(func=st.difference)

    # correlation(): returns a correlation matrix
    parser_stat_correlation = subparsers_stat.add_parser("correlation", help="Compute correlation matrix")

    parser_stat_correlation.add_argument("--df", type=str, help="Input file path",required=True)

    parser_stat_correlation.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_stat_correlation.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_correlation.csv')

    parser_stat_correlation.add_argument("--var_cols", nargs="+", help="List of 2 variable columns for tidy format")
    parser_stat_correlation.add_argument("--value_cols", nargs="+", help="List of numerical columns to correlate")
    parser_stat_correlation.add_argument("--method", type=str, default="pearson", choices=["pearson", "spearman", "kendall"],
                                         help="Correlation method to use (Default: pearson)")
    parser_stat_correlation.add_argument("--numeric_only", action="store_true", help="Only use numeric columns (Default: True)")

    parser_stat_correlation.set_defaults(func=st.correlation)

    # compare(): computes FC, pval, and log transformations relative to a specified condition
    parser_stat_compare = subparsers_stat.add_parser("compare", help="Compare conditions using FC, p-values, and log transforms")

    parser_stat_compare.add_argument("--df", type=str, help="Input file path",required=True)
    parser_stat_compare.add_argument("--sample", type=str, help="Sample column name",required=True)
    parser_stat_compare.add_argument("--cond", type=str, help="Condition column name",required=True)
    parser_stat_compare.add_argument("--cond_comp", type=str, help="Condition to compare against",required=True)
    parser_stat_compare.add_argument("--var", type=str, help="Variable column name",required=True)
    parser_stat_compare.add_argument("--count", type=str, help="Count column name",required=True)

    parser_stat_compare.add_argument("--psuedocount", type=int, default=1, help="Pseudocount to avoid log(0) or divide-by-zero errors")
    parser_stat_compare.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_stat_compare.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_compare.csv')

    parser_stat_compare.set_defaults(func=st.compare)

    '''
    edms.gen.io:
    - in_subs(): moves all files with a given suffix into subfolders named after the files (excluding the suffix).
    - out_subs(): recursively moves all files from subdirectories into the parent directory and delete the emptied subdirectories.
    - create_sh(): creates a shell script with SLURM job submission parameters for Harvard FASRC cluster.
    - split_R1_R2(): split paired reads into new R1 and R2 subdirectories at the parent directory
    - excel_csvs(): exports excel file to .csv files in specified directory  
    '''
    parser_io = subparsers.add_parser("io", help="Input/Output")
    subparsers_io = parser_io.add_subparsers()
    
    # Create subparsers for commands
    parser_io_in_subs = subparsers_io.add_parser("in_subs", help="*No FASRC* Moves all files with a given suffix into subfolders named after the files (excluding the suffix)")
    parser_io_out_subs = subparsers_io.add_parser("out_subs", help="*No FASRC* Delete subdirectories and move their files to the parent directory")
    parser_io_create_sh = subparsers_io.add_parser("create_sh", help='Generate SLURM shell script for Harvard FASRC cluster.')
    parser_io_split_R1_R2 = subparsers_io.add_parser("split_R1_R2", help='*No FASRC* Split paired reads into new R1 and R2 subdirectories at the parent directory.')
    parser_io_excel_csvs = subparsers_io.add_parser("excel_csvs", help='Exports excel file to .csv files in specified directory.')
    
    # Add common arguments
    for parser_io_common in [parser_io_in_subs,parser_io_out_subs,parser_io_split_R1_R2]:
        parser_io_common.add_argument("--dir", help="Path to parent directory", type=str, default='.')
    
    # in_subs() arguments
    parser_io_in_subs.add_argument("--suf", help="File suffix (e.g., '.txt', '.csv') to filter files.", type=int, required=True) 
    
    # create_sh(): creates a shell script with SLURM job submission parameters for Harvard FASRC cluster.
    parser_io_create_sh.add_argument('--dir', type=str, help='Directory to save the shell script.', default='.')
    parser_io_create_sh.add_argument('--file', type=str, help='Name of the shell script file to create.',default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.sh')
    parser_io_create_sh.add_argument('--cores', type=int, default=1, help='Number of CPU cores to request.')
    parser_io_create_sh.add_argument('--partition', type=str, default='serial_requeue', help='SLURM partition to use.')
    parser_io_create_sh.add_argument('--time', type=str, default='0-00:10', help='Job run time in D-HH:MM format.')
    parser_io_create_sh.add_argument('--mem', type=int, default=1000, help='Memory in MB.')
    parser_io_create_sh.add_argument('--email', type=str, default=None, help='Notification email address.')
    parser_io_create_sh.add_argument('--python', type=str, default='python/3.12.5-fasrc01', help='Python module to load.')
    parser_io_create_sh.add_argument('--env', type=str, default='edms', help='Conda environment to activate.')

    # excel_csvs(): exports excel file to .csv files in specified directory 
    parser_io_excel_csvs.add_argument('--pt', type=str, help='Excel file path', required=True)
    parser_io_excel_csvs.add_argument('--dir', type=str, help='Output directory path (Default: same directory as excel file name).',default='')

    # Call command functions
    parser_io_in_subs.set_defaults(func=io.in_subs)
    parser_io_out_subs.set_defaults(func=io.out_subs)
    parser_io_create_sh.set_defaults(func=io.create_sh)
    parser_io_split_R1_R2.set_defaults(func=io.split_R1_R2)
    parser_io_excel_csvs.set_defaults(func=io.excel_csvs)

    '''
    edms.gen.cli:
    - access(): make all files and subdirectories accessible on Harvard FASRC
    - smaller_fastq(): create new subdirectory containing fastqs with the # of reads limited
    - create_export_var(): create a persistent environment variable by adding it to the user's shell config.
    '''
    parser_cli = subparsers.add_parser("cli", help="Command Line Interaction")
    subparsers_cli = parser_cli.add_subparsers()
    
    # Create subparsers for commands
    parser_cli_access = subparsers_cli.add_parser("access", help="Make all files and subdirectories accessible on Harvard FASRC")
    parser_cli_smaller_fastq = subparsers_cli.add_parser("smaller_fastq", help="Create new subdirectory containing fastqs with the # of reads limited")
    parser_cli_create_export_var = subparsers_cli.add_parser("create_export_var", help="Create a persistent export variable by adding it to the user's shell config.")
    parser_cli_view_export_vars = subparsers_cli.add_parser("view_export_vars", help="View the current export variables in the user's shell config.")
    
    # Add common arguments
    for parser_cli_common in [parser_cli_access, parser_cli_smaller_fastq]:
        parser_cli_common.add_argument("--pt", help="Path to parent directory", type=str, default='.')
    
    # Smaller_fastq arguments
    parser_cli_smaller_fastq.add_argument("--reads", help="# of reads per fastq file", type=int, default='100000') 
    parser_cli_smaller_fastq.add_argument("--suf", help="Fastq file suffix", type=int, default=".fastq.gz") 
    
    # create_export_var arguments
    parser_cli_create_export_var.add_argument("--name", help="Name of the environment variable (e.g., MYPROJ)", required=True)
    parser_cli_create_export_var.add_argument("--pt", help="Path the variable should point to (e.g., ~/projects/myproj)", required=True)
    parser_cli_create_export_var.add_argument("--shell", choices=["bash", "zsh"], default=argparse.SUPPRESS, help="Shell type)")
    
    # view_export_var arguments
    parser_cli_view_export_vars.add_argument("--shell", choices=["bash", "zsh"], default=argparse.SUPPRESS, help="Shell type")

    # set default functions
    parser_cli_access.set_defaults(func=cli.access)
    parser_cli_smaller_fastq.set_defaults(func=cli.smaller_fastq)
    parser_cli_create_export_var.set_defaults(func=cli.create_export_var)
    parser_cli_view_export_vars.set_defaults(func=cli.view_export_vars)

    '''
    edms.dat.cosmic:
    - mutations(): returns COSMIC mutations dataframe for a given gene
    - cds_group(): plot COSMIC mutations histogram with CDS regions highlighted in different colors
    - priority_muts(): returns the shared sequences library dataframe with priority mutations
    - priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    - editor_mutations(): returns and plots editor accessible COSMIC mutations
    '''
    parser_cosmic = subparsers.add_parser("cosmic", help="COSMIC Database")
    subparsers_cosmic = parser_cosmic.add_subparsers()

    # mutations(): returns COSMIC mutations dataframe for a given gene
    parser_cosmic_mutations = subparsers_cosmic.add_parser("mutations", help="Extract COSMIC mutations")

    parser_cosmic_mutations.add_argument("--df", type=str, help="Input file path", required=True)

    parser_cosmic_mutations.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cosmic_mutations.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_cosmic_mutations.csv')

    parser_cosmic_mutations.set_defaults(func=co.mutations)
    
    # cds_group(): plot COSMIC mutations histogram with CDS regions highlighted in different colors
    parser_cds_group = subparsers_cosmic.add_parser("cds_group", help="Plot COSMIC mutation histogram with CDS regions highlighted")

    parser_cds_group.add_argument("--df_cosmic", type=str, help="COSMIC mutations() file path", required=True)
    parser_cds_group.add_argument("--df_cds", type=str, help="CDS region file path (with columns: gene, CDS, start, end)", required=True)

    parser_cds_group.add_argument("--out_dir", type=str, help="Output directory for plot",default='../out')

    parser_cds_group.set_defaults(func=co.cds_group)

    # priority_muts: returns the shared sequences library dataframe with priority mutations
    parser_cosmic_priority_muts = subparsers_cosmic.add_parser("priority_muts", help="Identify priority mutations in shared pegRNA library")

    parser_cosmic_priority_muts.add_argument("--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path", required=True)
    parser_cosmic_priority_muts.add_argument("--df_cosmic", type=str, help="COSMIC mutations() dataframe file path",required=True)

    parser_cosmic_priority_muts.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cosmic_priority_muts.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_shared_mutations.csv')

    parser_cosmic_priority_muts.set_defaults(func=co.priority_muts)

    # priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    parser_cosmic_priority_edits = subparsers_cosmic.add_parser("priority_edits", help="Identify clinically-relevant prime edits from shared pegRNA sequences")

    parser_cosmic_priority_edits.add_argument("--pegRNAs", type=str, help="pegRNAs library dataframe file path",required=True)
    parser_cosmic_priority_edits.add_argument("--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path",required=True)
    parser_cosmic_priority_edits.add_argument("--df_cosmic", type=str, help="COSMIC mutations() dataframe file path",required=True)
    
    parser_cosmic_priority_edits.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cosmic_priority_edits.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_priority.csv')

    parser_cosmic_priority_edits.set_defaults(func=co.priority_edits)

    # editor_mutations(): returns and plots editor accessible COSMIC mutations
    parser_editor_muts = subparsers_cosmic.add_parser("editor_mutations", help="Plot editor-accessible COSMIC mutations using BESCAN library")

    parser_editor_muts.add_argument("--df_cosmic", type=str, help="COSMIC mutations() dataframe file path",required=True)
    parser_editor_muts.add_argument("--df_bescan", type=str, help="BESCAN sgRNA library dataframe file path",required=True)

    parser_editor_muts.add_argument("--out_dir", type=str, help="Output directory for plots or results",default='../out')

    parser_editor_muts.set_defaults(func=co.editor_mutations)

    '''
    edms.dat.cvar:
    - mutations(): returns ClinVar mutations dataframe for a given gene
    - priority_muts: returns the shared sequences library dataframe with priority mutations
    - priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    '''
    parser_cvar = subparsers.add_parser("cvar", help="ClinVar Database")
    subparsers_cvar = parser_cvar.add_subparsers()

    # mutations(): returns ClinVar mutations dataframe for a given gene
    parser_cvar_mutations = subparsers_cvar.add_parser("mutations", help="Extract ClinVar mutations")

    parser_cvar_mutations.add_argument("--df", type=str, help="Input file path", required=True)
    parser_cvar_mutations.add_argument("--gene_name", type=str, help="Gene name", required=True)

    parser_cvar_mutations.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cvar_mutations.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_cvar_mutations.csv')

    parser_cvar_mutations.set_defaults(func=cvar.mutations)

    # priority_muts: returns the shared sequences library dataframe with priority mutations
    parser_cvar_priority_muts = subparsers_cvar.add_parser("priority_muts", help="Identify priority mutations in shared pegRNA library")

    parser_cvar_priority_muts.add_argument("--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path", required=True)
    parser_cvar_priority_muts.add_argument("--df_clinvar", type=str, help="ClinVar mutations() dataframe file path",required=True)

    parser_cvar_priority_muts.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cvar_priority_muts.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_shared_mutations.csv')

    parser_cvar_priority_muts.set_defaults(func=cvar.priority_muts)

    # priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    parser_cvar_priority_edits = subparsers_cvar.add_parser("priority_edits", help="Identify clinically-relevant prime edits from shared pegRNA sequences")

    parser_cvar_priority_edits.add_argument("--pegRNAs", type=str, help="pegRNAs library dataframe file path",required=True)
    parser_cvar_priority_edits.add_argument("--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path",required=True)
    parser_cvar_priority_edits.add_argument("--df_clinvar", type=str, help="ClinVar mutations() dataframe file path",required=True)

    parser_cvar_priority_edits.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cvar_priority_edits.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_priority.csv')

    parser_cvar_priority_edits.set_defaults(func=cvar.priority_edits)

    '''
    edms.bio.ngs:
    - pcrs(): generates NGS PCR plan automatically
    - umis(): determine the ug of gDNA and reads required for genotyping with UMIs.
    - compute_distance_matrix(): compute pairwise Hamming distance matrix for a list of sequences stored in a dataframe
    '''
    parser_ngs = subparsers.add_parser("ngs", help="Next generation sequencing")
    subparsers_ngs = parser_ngs.add_subparsers()

    # pcrs(): generates NGS PCR plan automatically
    parser_ngs_pcrs = subparsers_ngs.add_parser("pcrs", help="Plan NGS PCRs")
    
    # pcrs(): Core parameters
    parser_ngs_pcrs.add_argument("--df", help="Input file", type=str, required=True)
    
    parser_ngs_pcrs.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_ngs_pcrs.add_argument("--file", help="Output file name (.xlsx)", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_NGS_plan.xlsx')
    parser_ngs_pcrs.add_argument("--ultra", help="Using NEB Ultra II reagents", action="store_true")
    parser_ngs_pcrs.add_argument("--pcr1_cycles", help="Number of cycles for PCR1", type=str, default='30')
    parser_ngs_pcrs.add_argument("--pcr2_cycles", help="Number of cycles for PCR2", type=str, default='8')
    parser_ngs_pcrs.add_argument("--umi_cycles", help="Number of cycles for PCR1", type=str, default='3')
    parser_ngs_pcrs.add_argument("--pcr1.5_Tm",dest='pcr1_5_Tm', help="Annealing temperature for PCR1.5", type=str, default='65')
    parser_ngs_pcrs.add_argument('--pcr1_total_uL', type=int, default=20, help='PCR1 Total reaction volume (uL)')
    parser_ngs_pcrs.add_argument('--pcr2_total_uL', type=int, default=20, help='PCR2 Total reaction volume (uL)')
    parser_ngs_pcrs.add_argument('--mm_x', type=float, default=1.1, help='Master mix multiplier')
    
    # pcrs(): Column names
    parser_ngs_pcrs.add_argument('--gDNA_id_col', default='ID', help='gDNA ID column name')
    parser_ngs_pcrs.add_argument('--pcr1_id_col', default='PCR1 ID', help='PCR1 ID column name')
    parser_ngs_pcrs.add_argument('--pcr1_fwd_col', default='PCR1 FWD', help='PCR1 FWD column name')
    parser_ngs_pcrs.add_argument('--pcr1_rev_col', default='PCR1 REV', help='PCR1 REV column name')
    parser_ngs_pcrs.add_argument('--pcr2_id_col', default='PCR2 ID', help='PCR2 ID column name')
    parser_ngs_pcrs.add_argument('--pcr2_fwd_col', default='PCR2 FWD', help='PCR2 FWD column name')
    parser_ngs_pcrs.add_argument('--pcr2_rev_col', default='PCR2 REV', help='PCR2 REV column name')
    parser_ngs_pcrs.add_argument('--umi_col', default='UMI', help='UMI column name')

    # pcrs(): Stock concentrations
    parser_ngs_pcrs.add_argument('--Q5_mm_x_stock', type=float, default=5, help='Q5 reaction master mix stock (X)')
    parser_ngs_pcrs.add_argument('--dNTP_mM_stock', type=float, default=10, help='dNTP stock concentration (mM)')
    parser_ngs_pcrs.add_argument('--fwd_uM_stock', type=float, default=10, help='Forward primer stock concentration (uM)')
    parser_ngs_pcrs.add_argument('--rev_uM_stock', type=float, default=10, help='Reverse primer stock concentration (uM)')
    parser_ngs_pcrs.add_argument('--Q5_U_uL_stock', type=float, default=2, help='Q5 Polymerase stock (U/uL)')

    # pcrs(): Desired concentrations
    parser_ngs_pcrs.add_argument('--Q5_mm_x_desired', type=float, default=1, help='Q5 reaction master mix desired (X)')
    parser_ngs_pcrs.add_argument('--dNTP_mM_desired', type=float, default=0.2, help='dNTP desired concentration (mM)')
    parser_ngs_pcrs.add_argument('--fwd_uM_desired', type=float, default=0.5, help='Forward primer desired concentration (uM)')
    parser_ngs_pcrs.add_argument('--rev_uM_desired', type=float, default=0.5, help='Reverse primer desired concentration (uM)')
    parser_ngs_pcrs.add_argument('--Q5_U_uL_desired', type=float, default=0.02, help='Q5 Polymerase desired amount (U/uL)')

    parser_ngs_pcrs.set_defaults(func=ngs.pcrs)
    
    # umis(): determine the ug of gDNA and reads required for genotyping with UMIs.
    parser_ngs_umis = subparsers_ngs.add_parser("umis", help="Determine ug of gDNA and reads required for genotyping with UMIs")

    parser_ngs_umis.add_argument("--genotypes", help="# of intended genotypes per sample", type=int, required=True)

    parser_ngs_umis.add_argument("--samples", help="# of samples to be processed (Default: 1", type=int, default=1)
    parser_ngs_umis.add_argument("--cell_coverage", help="Desired coverage per genotype (Default: 1000)", type=int, default=1000)
    parser_ngs_umis.add_argument("--ug_gDNA_per_cell", help="Amount of genomic DNA per cell in micrograms (Default: 6x10^(-6) ug/cell)", type=float, default=6*10**-6)
    parser_ngs_umis.add_argument("--ploidy_per_cell", help="Ploidy level of the cells (Default: 2 for diploid)", type=int, default=2)
    parser_ngs_umis.add_argument("--umi_coverage", help="Average # of reads per UMI (Default: 5)", type=int, default=5)

    parser_ngs_umis.set_defaults(func=ngs.umis)
    
    # hamming_distance_matrix(): compute pairwise Hamming distance matrix for a list of sequences stored in a dataframe
    parser_ngs_hamming = subparsers_ngs.add_parser("hamming", help="Compute pairwise Hamming distance matrix")
    
    parser_ngs_hamming.add_argument("--df", help="Input file", type=str, required=True)
    parser_ngs_hamming.add_argument("--id", help="ID column name", type=str, required=True)
    parser_ngs_hamming.add_argument("--seqs", help="Sequences column name", type=str, required=True)
    
    parser_ngs_hamming.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_ngs_hamming.add_argument("--file", help="Output file name", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_hamming.csv')
    
    parser_ngs_hamming.set_defaults(func=ngs.hamming_distance_matrix)

    '''
    edms.bio.sanger:
    - pcrs(): generates Sanger PCR plan automatically
    '''
    parser_sanger = subparsers.add_parser("sanger", help="Sanger sequencing")
    subparsers_sanger = parser_sanger.add_subparsers()

    # pcrs(): generates Sanger PCR plan automatically
    parser_sanger_pcrs = subparsers_sanger.add_parser("pcrs", help="Plan Sanger PCRs")
    
    # pcrs(): Core parameters
    parser_sanger_pcrs.add_argument("--df", help="Input file", type=str, required=True)
    parser_sanger_pcrs.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_sanger_pcrs.add_argument("--file", help="Output file name (.xlsx)", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_Sanger_plan.xlsx')
    parser_sanger_pcrs.add_argument("--cycles", help="Number of cycles for PCR1", type=str, default='30')
    parser_sanger_pcrs.add_argument("--ultra", help="Using NEB Ultra II reagents", action="store_true")
    parser_sanger_pcrs.add_argument('--total_uL', type=int, default=20, help='Total reaction volume (uL)')
    parser_sanger_pcrs.add_argument('--mm_x', type=float, default=1.1, help='Master mix multiplier')
    
    # pcrs(): Column names
    parser_sanger_pcrs.add_argument('--gDNA_id_col', default='ID', help='gDNA ID column name')
    parser_sanger_pcrs.add_argument('--pcr1_id_col', default='PCR1 ID', help='PCR1 ID column name')
    parser_sanger_pcrs.add_argument('--pcr1_fwd_col', default='PCR1 FWD', help='PCR1 FWD column name')
    parser_sanger_pcrs.add_argument('--pcr1_rev_col', default='PCR1 REV', help='PCR1 REV column name')
   
    # pcrs(): Stock concentrations
    parser_sanger_pcrs.add_argument('--Q5_mm_x_stock', type=float, default=5, help='Q5 reaction master mix stock (X)')
    parser_sanger_pcrs.add_argument('--dNTP_mM_stock', type=float, default=10, help='dNTP stock concentration (mM)')
    parser_sanger_pcrs.add_argument('--fwd_uM_stock', type=float, default=10, help='Forward primer stock concentration (uM)')
    parser_sanger_pcrs.add_argument('--rev_uM_stock', type=float, default=10, help='Reverse primer stock concentration (uM)')
    parser_sanger_pcrs.add_argument('--Q5_U_uL_stock', type=float, default=2, help='Q5 Polymerase stock (U/uL)')

    # pcrs(): Desired concentrations
    parser_sanger_pcrs.add_argument('--Q5_mm_x_desired', type=float, default=1, help='Q5 reaction master mix desired (X)')
    parser_sanger_pcrs.add_argument('--dNTP_mM_desired', type=float, default=0.2, help='dNTP desired concentration (mM)')
    parser_sanger_pcrs.add_argument('--fwd_uM_desired', type=float, default=0.5, help='Forward primer desired concentration (uM)')
    parser_sanger_pcrs.add_argument('--rev_uM_desired', type=float, default=0.5, help='Reverse primer desired concentration (uM)')
    parser_sanger_pcrs.add_argument('--Q5_U_uL_desired', type=float, default=0.02, help='Q5 Polymerase desired amount (U/uL)')

    parser_sanger_pcrs.set_defaults(func=sanger.pcrs)
    
    '''
    edms.bio.clone
    - sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
    - epegRNAs(): design GG cloning oligonucleotides for prime editing epegRNAs
    - ngRNAs(): design GG cloning oligonucleotides for prime editing ngRNAs
    - epegRNA_pool(): makes twist oligonucleotides for prime editing
    - pcr_sim(): returns dataframe with simulated pcr product
    - off_targets(): Find off-target sequences for a list of sequences using pairwise alignment.
    '''
    parser_clone = subparsers.add_parser("clone", help="Molecular cloning")
    subparsers_clone = parser_clone.add_subparsers()

    # sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
    parser_clone_sgRNAs = subparsers_clone.add_parser("sgRNAs", help="Design GG oligos for sgRNAs (cutting or BE)")
    
    parser_clone_sgRNAs.add_argument("--df", type=str, help="Input file path",required=True)
    parser_clone_sgRNAs.add_argument("--id", type=str, help="Column name for unique sgRNA identifier",required=True)

    parser_clone_sgRNAs.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_clone_sgRNAs.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_sgRNAs.csv')

    parser_clone_sgRNAs.add_argument("--dont_tG", dest="tG", default=True, action="store_false", help="Don't add 5' G to spacer if needed")
    parser_clone_sgRNAs.add_argument("--dont_order", dest="order", default=True, action="store_false", help="Don't format output for ordering oligos")
    parser_clone_sgRNAs.add_argument("--spacer", type=str, default="Spacer_sequence", help="Column name for spacer sequence")
    parser_clone_sgRNAs.add_argument("--t5", type=str, default="CACC", help="Top oligo 5' overhang")
    parser_clone_sgRNAs.add_argument("--t3", type=str, default="", help="Top oligo 3' overhang")
    parser_clone_sgRNAs.add_argument("--b5", type=str, default="AAAC", help="Bottom oligo 5' overhang (revcom)")
    parser_clone_sgRNAs.add_argument("--b3", type=str, default="", help="Bottom oligo 3' overhang (revcom)")

    parser_clone_sgRNAs.set_defaults(func=cl.sgRNAs)

    # epegRNAs(): design GG cloning oligonucleotides for prime editing epegRNAs
    parser_clone_epegRNAs = subparsers_clone.add_parser("epegRNAs", help="Design GG oligos for epegRNAs")

    parser_clone_epegRNAs.add_argument("--df", type=str, help="Input file path", required=True)
    parser_clone_epegRNAs.add_argument("--id", type=str, help="Column name for unique sequence identifier",required=True)

    parser_clone_epegRNAs.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_clone_epegRNAs.add_argument("--file", help="Output file name", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_epegRNAs.csv')

    parser_clone_epegRNAs.add_argument("--dont_tG", dest="tG", default=True, action="store_false", help="Don't add 5' G to spacer if needed")
    parser_clone_epegRNAs.add_argument("--dont_order", dest="order", default=True, action="store_false", help="Don't format output for ordering oligos")
    parser_clone_epegRNAs.add_argument("--order_scaffold", action="store_true", help="Include scaffold sequence in the oligo order")
    parser_clone_epegRNAs.add_argument("--dont_make_extension", dest="make_extension", default=True, action="store_false", help="Don't build extension from RTT, PBS, and linker")
    parser_clone_epegRNAs.add_argument("--spacer", type=str, default="Spacer_sequence", help="Column name for spacer sequence")
    parser_clone_epegRNAs.add_argument("--spacer_t5", type=str, default="CACC", help="Top 5' overhang for spacer")
    parser_clone_epegRNAs.add_argument("--spacer_t3", type=str, default="GTTTAAGAGC", help="Top 3' overhang for spacer")
    parser_clone_epegRNAs.add_argument("--spacer_b5", type=str, default="", help="Bottom 5' overhang for spacer")
    parser_clone_epegRNAs.add_argument("--spacer_b3", type=str, default="", help="Bottom 3' overhang for spacer")
    parser_clone_epegRNAs.add_argument("--extension", type=str, default="Extension_sequence", help="Column name for extension sequence")
    parser_clone_epegRNAs.add_argument("--extension_t5", type=str, default="", help="Top 5' overhang for extension")
    parser_clone_epegRNAs.add_argument("--extension_t3", type=str, default="", help="Top 3' overhang for extension")
    parser_clone_epegRNAs.add_argument("--extension_b5", type=str, default="CGCG", help="Bottom 5' overhang for extension")
    parser_clone_epegRNAs.add_argument("--extension_b3", type=str, default="GCACCGACTC", help="Bottom 3' overhang for extension")
    parser_clone_epegRNAs.add_argument("--RTT", type=str, default="RTT_sequence", help="Column name for RTT (reverse transcriptase template)")
    parser_clone_epegRNAs.add_argument("--PBS", type=str, default="PBS_sequence", help="Column name for PBS (primer binding site)")
    parser_clone_epegRNAs.add_argument("--linker", type=str, default="Linker_sequence", help="Column name for linker")

    parser_clone_epegRNAs.set_defaults(func=cl.epegRNAs)
    
    # ngRNAs(): design GG cloning oligonucleotides for prime editing ngRNAs
    parser_clone_ngRNAs = subparsers_clone.add_parser("ngRNAs", help="Design GG oligos for ngRNAs")

    parser_clone_ngRNAs.add_argument("--df", type=str, help="Input file path", required=True)
    parser_clone_ngRNAs.add_argument("--id", type=str, help="Column name for unique sequence identifier",required=True)

    parser_clone_ngRNAs.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_clone_ngRNAs.add_argument("--file", help="Output file name", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_epegRNAs.csv')
    
    parser_clone_ngRNAs.add_argument("--dont_tG", dest="tG", default=True, action="store_false", help="Don't add 5' G to spacer if needed")
    parser_clone_ngRNAs.add_argument("--dont_order", dest="order", default=True, action="store_false", help="Don't format output for ordering oligos")
    parser_clone_ngRNAs.add_argument("--order_scaffold", action="store_true", help="Include scaffold sequence in the oligo order")
    parser_clone_ngRNAs.add_argument("--spacer", type=str, default="Spacer_sequence", help="Column name for spacer sequence")
    parser_clone_ngRNAs.add_argument("--spacer_t5", type=str, default="CACC", help="Top strand 5' overhang")
    parser_clone_ngRNAs.add_argument("--spacer_t3", type=str, default="GTTTAAGAGC", help="Top strand 3' overhang")
    parser_clone_ngRNAs.add_argument("--spacer_b5", type=str, default="", help="Bottom strand 5' overhang")
    parser_clone_ngRNAs.add_argument("--spacer_b3", type=str, default="", help="Bottom strand 3' overhang")

    parser_clone_ngRNAs.set_defaults(func=cl.ngRNAs)

    # epegRNA_pool(): makes twist oligonucleotides for prime editing
    parser_clone_epegRNA_pool = subparsers_clone.add_parser("epegRNA_pool", help="Design GG oligos for pooled epegRNAs")
    
    parser_clone_epegRNA_pool.add_argument("--df", type=str, help="Input file path", required=True)

    parser_clone_epegRNA_pool.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_clone_epegRNA_pool.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_epegRNA_pool.csv')

    parser_clone_epegRNA_pool.add_argument("--tG", action="store_true", help="Add 5' G to spacers if needed")
    parser_clone_epegRNA_pool.add_argument("--make_extension", action="store_true", help="Build extension from RTT, PBS, and linker")
    parser_clone_epegRNA_pool.add_argument("--UMI_df", type=str, help="UMI sequences file path", default=argparse.SUPPRESS)
    parser_clone_epegRNA_pool.add_argument("--PCR_df", type=str, help="PCR primer and subpool barcode file path", default=argparse.SUPPRESS)
    parser_clone_epegRNA_pool.add_argument("--RE_type_IIS_df", type=str, help="RE Type IIS file path", default=argparse.SUPPRESS)
    parser_clone_epegRNA_pool.add_argument("--UMI_i", type=int, help="UMI start index (Default: 0)", default=0)
    parser_clone_epegRNA_pool.add_argument("--enzymes", type=str, nargs="+", help="List of Type IIS restriction enzymes to check for (Default: Esp3I)", default='Esp3I')
    parser_clone_epegRNA_pool.add_argument("--barcode", type=str, help="subpool barcode column name (Default: Barcode)", default='Barcode')
    parser_clone_epegRNA_pool.add_argument("--barcode_i", type=int, help="subpool barcode start index (Default: 0)", default=0)
    parser_clone_epegRNA_pool.add_argument("--fwd_barcode_t5", type=str, default="Forward Barcode", help="Forward barcode column name")
    parser_clone_epegRNA_pool.add_argument("--rev_barcode_t3", type=str, default="Reverse Barcode", help="Reverse barcode column name")
    parser_clone_epegRNA_pool.add_argument("--Esp3I_hU6", type=str, default="Esp3I_hU6", help="Esp3I_hU6 column name")
    parser_clone_epegRNA_pool.add_argument("--tevopreQ1_Esp3I", type=str, default="tevopreQ1_Esp3I", help="tevopreQ1_Esp3I column name")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_spacer", type=str, default="Spacer_sequence", help="epegRNA spacer column")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_scaffold", type=str, default="Scaffold_sequence", help="epegRNA scaffold column")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_extension", type=str, default="Extension_sequence", help="epegRNA extension column")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_RTT", type=str, default="RTT_sequence", help="epegRNA RTT column name")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_PBS", type=str, default="PBS_sequence", help="epegRNA PBS column name")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_linker", type=str, default="Linker_sequence", help="epegRNA Linker column name")

    parser_clone_epegRNA_pool.set_defaults(func=cl.epegRNA_pool) 
    
    # umi(): generates unique molecular identifiers (UMIs) of specified length, GC content, and Hamming distance
    parser_clone_umi = subparsers_clone.add_parser("umi", help="Generate unique molecular identifiers (UMIs) of specified length, GC content, and Hamming distance")

    parser_clone_umi.add_argument("--length", type=int, help="Length of UMI (Default: 15)", default=15)
    parser_clone_umi.add_argument("--GC_fract", type=parse_tuple_float, default=(0.4,0.6), help="Pair of GC content boundaries written as fractions (Default: 0.4,0.6)")
    parser_clone_umi.add_argument("--hamming", type=int, help="Minimum Hamming distance between UMIs (Default: 4)", default=4)
    parser_clone_umi.add_argument("--nrows", type=int, help="# of UMIs to compare iteratively for hamming filtering (Default: 1000)", default=1000)
    parser_clone_umi.add_argument("--pt", type=str, help="Shuffled UMI file path if already made (Default: None)", default=argparse.SUPPRESS)
    parser_clone_umi.add_argument("--dir", type=str, help="Output directory", default='../out')

    parser_clone_umi.set_defaults(func=cl.umi)

    # pcr_sim(): returns dataframe with simulated pcr product 
    parser_clone_pcrsim = subparsers_clone.add_parser("pcr_sim", help="Simulate PCR product from template and primer sequences")

    parser_clone_pcrsim.add_argument("--df", type=str, help="Input dataframe or file path containing template and primers", required=True)
    parser_clone_pcrsim.add_argument("--template_col", type=str, help="Column name for template sequence", required=True)
    parser_clone_pcrsim.add_argument("--fwd_bind_col", type=str, help="Column name for forward primer binding region", required=True)
    parser_clone_pcrsim.add_argument("--rev_bind_col", type=str, help="Column name for reverse primer binding region", required=True)

    parser_clone_pcrsim.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_clone_pcrsim.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pcr_sim.csv')

    parser_clone_pcrsim.add_argument("--fwd_ext_col", type=str, help="Column name for forward primer extension region")
    parser_clone_pcrsim.add_argument("--rev_ext_col", type=str, help="Column name for reverse primer extension region")
    parser_clone_pcrsim.add_argument("--product_col", type=str, default="PCR Product", help="Column name for output PCR product")
    
    parser_clone_pcrsim.set_defaults(func=cl.pcr_sim)

    # off_targets(): Find off-target sequences for a list of sequences using pairwise alignment.
    parser_clone_off_targets = subparsers_clone.add_parser("off_targets", help="Find off-target sequences for a list of sequences of the same length using pairwise alignment")

    parser_clone_off_targets.add_argument("--df", type=str, help="Input file path containing sequences", required=True)
    parser_clone_off_targets.add_argument("--col", type=str, help="Column name containing sequences to align", required=True)

    parser_clone_off_targets.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_clone_off_targets.add_argument("--ckpt", type=int, help="Checkpoint interval for saving progress (Default: 100)", default=100)
    parser_clone_off_targets.add_argument("--match_score", type=float, help="Match score for pairwise alignment", default=argparse.SUPPRESS)
    parser_clone_off_targets.add_argument("--mismatch_score", type=float, help="Mismatch score for pairwise alignment", default=argparse.SUPPRESS)
    parser_clone_off_targets.add_argument("--open_gap_score", type=float, help="Open gap score for pairwise alignment", default=argparse.SUPPRESS)
    parser_clone_off_targets.add_argument("--extend_gap_score", type=float, help="Extend gap score for pairwise alignment", default=argparse.SUPPRESS)

    parser_clone_off_targets.set_defaults(func=cl.off_targets)

    '''
    edms.bio.transfect
    - PE3(): generates PE3 transfection plan for HEK293T cells (Default: 96-well plate in triplicate using L2000)
    - virus(): generates transfection plan for virus production from HEK293T cells (Default: 6-well plate using L3000)
    '''
    parser_transfect = subparsers.add_parser("transfect", help="Transfection")
    subparsers_transfect = parser_transfect.add_subparsers()

    # PE3(): generates PE3 transfection plan for HEK293T cells (Default: 96-well plate in triplicate using L2000)
    parser_transfect_PE3 = subparsers_transfect.add_parser("PE3", help="Plan PE3 transfection")
    
    parser_transfect_PE3.add_argument("--plasmids", type=str, help="Path to plasmids file", required=True)
    parser_transfect_PE3.add_argument("--epegRNAs", type=str, help="Path to epegRNAs file", required=True)
    parser_transfect_PE3.add_argument("--ngRNAs", type=str, help="Path to ngRNAs file", required=True)

    parser_transfect_PE3.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_transfect_PE3.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_transfect_PE3.csv')

    parser_transfect_PE3.add_argument("--pegRNA_number_col", type=str, default="pegRNA_number", help="Column name for pegRNA number")
    parser_transfect_PE3.add_argument("--epegRNAs_name_col", type=str, default="Name", help="Column name for epegRNA name")
    parser_transfect_PE3.add_argument("--ngRNAs_name_col", type=str, default="Name", help="Column name for ngRNA name")
    parser_transfect_PE3.add_argument("--plasmid_col", type=str, default="Plasmid", help="Column name for plasmid name")
    parser_transfect_PE3.add_argument("--description_col", type=str, default="Description", help="Column name for plasmid description")
    parser_transfect_PE3.add_argument("--colony_col", type=str, default="Colony", help="Column name for colony name")
    parser_transfect_PE3.add_argument("--ng_uL_col", type=str, default="ng/uL", help="Column name for ng/uL concentration")
    parser_transfect_PE3.add_argument("--PE_plasmid", type=str, default="pMUZ86.7", help="Name of PE plasmid to search for")
    parser_transfect_PE3.add_argument("--reps", type=int, default=3, help="Number of replicates")
    parser_transfect_PE3.add_argument("--mm_x", type=float, default=1.1, help="Master mix multiplier")
    parser_transfect_PE3.add_argument("--epegRNA_ng", type=int, default=66, help="ng of epegRNA per well")
    parser_transfect_PE3.add_argument("--ngRNA_ng", type=int, default=22, help="ng of ngRNA per well")
    parser_transfect_PE3.add_argument("--PE_ng", type=int, default=200, help="ng of PE plasmid per well")
    parser_transfect_PE3.add_argument("--well_uL", type=int, default=10, help="Total uL per well")

    parser_transfect_PE3.set_defaults(func=tf.PE3)

    # virus(): generates transfection plan for virus production from HEK293T cells (Default: 6-well plate using L3000)
    parser_transfect_virus = subparsers_transfect.add_parser("virus", help="Plan virus transfection")

    parser_transfect_virus.add_argument("--plasmids", type=str, help="Path to plasmids file", required=True)

    parser_transfect_virus.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_transfect_virus.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_transfect_virus.csv')

    parser_transfect_virus.add_argument("--plasmid_col", type=str, default="Plasmid", help="Column name for plasmid name")
    parser_transfect_virus.add_argument("--description_col", type=str, default="Description", help="Column name for plasmid description")
    parser_transfect_virus.add_argument("--colony_col", type=str, default="Colony", help="Column name for colony")
    parser_transfect_virus.add_argument("--ng_uL_col", type=str, default="ng/uL", help="Column name for ng/uL concentration")
    parser_transfect_virus.add_argument("--VSVG_plasmid", type=str, default="pMUZ26.6", help="Name of VSVG plasmid")
    parser_transfect_virus.add_argument("--GagPol_plasmid", type=str, default="pMUZ26.7", help="Name of GagPol plasmid")
    parser_transfect_virus.add_argument("--reps", type=int, default=1, help="Number of replicates")
    parser_transfect_virus.add_argument("--mm_x", type=float, default=1.1, help="Master mix multiplier")
    parser_transfect_virus.add_argument("--VSVG_ng", type=int, default=750, help="VSVG ng per well")
    parser_transfect_virus.add_argument("--GagPol_ng", type=int, default=1500, help="GagPol ng per well")
    parser_transfect_virus.add_argument("--transfer_ng", type=int, default=750, help="Transfer plasmid ng per well")
    parser_transfect_virus.add_argument("--well_uL", type=int, default=500, help="Total uL per well")

    parser_transfect_virus.set_defaults(func=tf.virus)

    '''
    edms.bio.qPCR:
    - ddCq(): computes ΔΔCq mean and error for all samples holding target pairs constant
    '''
    # ddCq(): computes ΔΔCq mean and error for all samples holding target pairs constant
    parser_ddcq = subparsers.add_parser("ddCq", help="Compute ΔΔCq values for RT-qPCR data")
    
    parser_ddcq.add_argument("--data", type=str, help="Input Cq file from CFX instrument",required=True)

    parser_ddcq.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_ddcq.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_qPCR_ddCq.csv')

    parser_ddcq.add_argument("--sample_col", type=str, default="Sample", help="Column name for sample ID")
    parser_ddcq.add_argument("--target_col", type=str, default="Target", help="Column name for target gene ID")
    parser_ddcq.add_argument("--Cq_col", type=str, default="Cq", help="Column name for Cq values")

    parser_ddcq.set_defaults(func=qPCR.ddCq)

    '''
    edms.bio.fastq:
    - revcom_fastqs(): write reverse complement of fastqs to a new directory
    - unzip_fastqs(): Unzip gzipped fastqs and write to a new directory
    - comb_fastqs(): Combines one or more (un)compressed fastqs files into a single (un)compressed fastq file
    - genotyping(): quantify edit outcomes workflow
    - abundances(): quantify desired edits count & fraction per sample
    - count_motif(): returns a dataframe with the sequence motif location with mismatches per read for every fastq file in a directory
    - plot_motif(): generate plots highlighting motif mismatches, locations, and sequences
    - plot_alignments(): generate line & distribution plots from fastq alignments dictionary
    - count_region(): align read region from fastq directory to the annotated library with mismatches; plot and return fastq alignments dictionary
    - count_alignments(): align reads from fastq directory to annotated library with mismatches; plot and return fastq alignments dictionary
    - plot_paired(): generate stacked bar plots from paired_regions() dataframe
    - paired_regions(): quantify, plot, & return (un)paired regions that aligned to the annotated library
    - count_signatures(): generate signatures from fastq read region alignments to WT sequence; count signatures, plot and return fastq signatures dataframe
    - editing_per_library(): Determine editing relative library abundance
    - extract_umis(): extract UMIs using umi_tools
    - trim_motifs(): trimming motifs with cutadapt
    - make_sams(): generates alignments saved as a SAM files using bowtie2
    - make_bams(): converts SAM files to BAM files using samtools
    - bam_umi_tags(): copy UMI in read ID to RX tag in BAM files using fgbio
    - group_umis(): group BAM files by UMI using fgbio
    - consensus_umis(): generate consensus sequences from grouped UMIs using fgbio
    - bam_to_fastq(): convert BAM files to FASTQ files using samtools
    '''
    parser_fastq = subparsers.add_parser("fastq", help="FASTQ files")
    subparsers_fastq = parser_fastq.add_subparsers()

    parser_fastq_revcom = subparsers_fastq.add_parser("revcom", help="Reverse complement all FASTQ files in a directory")
    parser_fastq_unzip = subparsers_fastq.add_parser("unzip", help="Unzip gzipped FASTQ files to a new directory")
    parser_fastq_comb = subparsers_fastq.add_parser("comb", help="Combine multiple FASTQ files into a single FASTQ (.fastq or .fastq.gz)")
    parser_fastq_genotyping = subparsers_fastq.add_parser("genotyping", help="Quantify edit outcomes workflow")
    parser_fastq_abundances = subparsers_fastq.add_parser("abundances", help="Quantify edit outcomes count & fraction per sample")
    parser_fastq_count_motif = subparsers_fastq.add_parser("count_motif", help="Count motif occurrences in FASTQ files")
    parser_fastq_plot_motif = subparsers_fastq.add_parser("plot_motif", help="Plot motif occurrences from FASTQ files")
    parser_fastq_plot_alignments = subparsers_fastq.add_parser("plot_alignments", help="Plot alignments from FASTQ files")
    parser_fastq_count_region = subparsers_fastq.add_parser("count_region", help="Count region occurrences in FASTQ files")
    parser_fastq_count_alignments = subparsers_fastq.add_parser("count_alignments", help="Count alignments in FASTQ files")
    parser_fastq_plot_paired = subparsers_fastq.add_parser("plot_paired", help="Plot paired regions from FASTQ files")
    parser_fastq_paired_regions = subparsers_fastq.add_parser("paired_regions", help="Extract paired regions from FASTQ files")
    parser_fastq_count_signatures = subparsers_fastq.add_parser("count_signatures", help="Generate signatures from fastq read region alignments to WT sequence")
    parser_fastq_editing_per_library = subparsers_fastq.add_parser("editing_per_library", help="Determine editing relative library abundance")
    parser_fastq_extract_umis = subparsers_fastq.add_parser("extract_umis", help="Extract UMIs using umi_tools")
    parser_fastq_trim_motifs = subparsers_fastq.add_parser("trim_motifs", help="Trim motifs with cutadapt")
    parser_fastq_make_sams = subparsers_fastq.add_parser("make_sams", help="Generate alignments saved as SAM files using bowtie2")
    parser_fastq_make_bams = subparsers_fastq.add_parser("make_bams", help="Convert SAM files to BAM files using samtools")
    parser_fastq_bam_umi_tags = subparsers_fastq.add_parser("bam_umi_tags", help="Copy UMI in read ID to RX tag in BAM files using fgbio")
    parser_fastq_group_umis = subparsers_fastq.add_parser("group_umis", help="Group BAM files by UMI using fgbio")
    parser_fastq_consensus_umis = subparsers_fastq.add_parser("consensus_umis", help="Generate consensus sequences from grouped UMIs using fgbio")
    parser_fastq_bam_to_fastq = subparsers_fastq.add_parser("bam_to_fastq", help="Convert BAM files to FASTQ files using samtools")

    # Add common arguments: revcom_fastqs(), unzip_fastqs(), comb_fastqs(), and genotyping()
    for parser_fastq_common in [parser_fastq_revcom,parser_fastq_unzip,parser_fastq_comb,parser_fastq_genotyping]:
        parser_fastq_common.add_argument("--in_dir", type=str, help="Input directory containing FASTQ files",default='.')
        parser_fastq_common.add_argument("--out_dir", type=str, help="Output directory",default = f'../out')

    # Add specific arguments: comb_fastqs()
    parser_fastq_comb.add_argument("--out_file", type=str, help="Name of output FASTQ file (.fastq or .fastq.gz)", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_comb.fastq.gz')

    # Add specific arguments: genotyping()
    parser_fastq_genotyping.add_argument("--out_file_prefix", type=str, help="Name of output file prefix", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    
    parser_fastq_genotyping.add_argument("--config_key", type=str, help="Configuration file key for sequence [flank5(genotype region)flank3} & res (first amino acid in genotype region)", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--sequence", type=str, help="Formatted sequence: flank5(genotype region)flank3", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--res", type=int, help="First amino acid number in genotype region", default=argparse.SUPPRESS)
    
    parser_fastq_genotyping.add_argument("--qall", type=int, help="Minimum Phred quality score for all bases", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--qtrim", type=int, help="Phred quality threshold for end trimming", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--qavg", type=int, help="Minimum average Phred quality score", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--qmask", type=int, help="Phred quality threshold for masking to N", default=argparse.SUPPRESS)

    parser_fastq_genotyping.add_argument("--save", action="store_true", help="Save read statistics and genotypes files", dest="save", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--no_save", action="store_false", help="Don't save read statistics and genotypes files", dest="save", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--masks", action="store_true", help="Include masked sequence and translation",default=False)
    parser_fastq_genotyping.add_argument("--keepX", action="store_true", help="Keep unknown translation (X) in output", default=False)

    parser_fastq_genotyping.add_argument("--match_score", type=float, help="Match score for pairwise alignment", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--mismatch_score", type=float, help="Mismatch score for pairwise alignment", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--open_gap_score", type=float, help="Open gap score for pairwise alignment", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--extend_gap_score", type=float, help="Extend gap score for pairwise alignment", default=argparse.SUPPRESS)
    
    # abundances():
    parser_fastq_abundances.add_argument("--df", help="Input file with sample, edit, count, & fraction information", required=True)
    parser_fastq_abundances.add_argument("--desired_edits", nargs="+", help="List of desired edits to isolate (space-separated)",required=True)

    parser_fastq_abundances.add_argument("--sample_col", default="sample", help="Column for sample ID (Default: 'sample')")
    parser_fastq_abundances.add_argument("--edit_col", default="edit", help="Column for edit identifier (Default: 'edit')")
    parser_fastq_abundances.add_argument("--count_col", default="count", help="Column for edit count (Default: 'count')")
    parser_fastq_abundances.add_argument("--fraction_col", default="fraction", help="Column for edit fraction (Default: 'fraction')")
    parser_fastq_abundances.add_argument("--combinations", default=1, help="Maximum # of desired edit combinations to search for (Default: 1 => single edits)")

    # count_motif():
    parser_fastq_count_motif.add_argument("--fastq_dir", help="Path to directory containing FASTQ files", required=True)
    parser_fastq_count_motif.add_argument("--pattern", help="Motif sequence pattern to search for", required=True)
    parser_fastq_count_motif.add_argument("--out_dir", help="Output directory to save results", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')

    parser_fastq_count_motif.add_argument("--motif", default="motif", help="Name of the motif (Default: 'motif')")
    parser_fastq_count_motif.add_argument("--max_distance", type=int, default=0, help="Maximum Levenshtein distance allowed (i.e., # of mismatches, Default: 0)")
    parser_fastq_count_motif.add_argument("--max_reads", type=int, default=0, help="Maximum # of reads to process per file")
    parser_fastq_count_motif.add_argument("--meta", type=str, help="Optional path to metadata CSV/TSV file with 'fastq_file' column")

    # plot_motif():
    parser_fastq_plot_motif.add_argument("--df", help="Path to count_motif() output file", required=True)

    parser_fastq_plot_motif.add_argument("--out_dir", type=str, help="Directory to save plots", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_plot_motif.add_argument("--plot_suf", type=str, default=".pdf", help="Plot file suffix (e.g. .pdf, .png)")
    parser_fastq_plot_motif.add_argument("--numeric", choices=["count", "fraction"], default="count", help="Numeric column to use for plotting (Default: 'count')")
    parser_fastq_plot_motif.add_argument("--id_col", default="fastq_file", help="Column used for sample ID (Default: 'fastq_file')")
    parser_fastq_plot_motif.add_argument("--id_axis", default="fastq", help="Label to use on the plot axis (Default: 'fastq')")
    parser_fastq_plot_motif.add_argument("--stack_figsize", type=parse_tuple_int, default=(7, 3), help="Stacked plot figure size formatted as 'width,height'")
    parser_fastq_plot_motif.add_argument("--heat_figsize", type=parse_tuple_int, help="Heatmap figure size formateed as 'width,height'")
    parser_fastq_plot_motif.add_argument("--cutoff_frac", type=float, default=0.01, help="Minimum fraction to include in y-axis (Default: 0.01)")

    # plot_alignments():
    parser_fastq_plot_alignments.add_argument("--fastq_alignments", help="Directory with the fastq alignments dictionary", required=True)
    parser_fastq_plot_alignments.add_argument("--align_col", help="Align column name in the annotated library reference file", required=True)
    parser_fastq_plot_alignments.add_argument("--id_col", help="ID column name in the annotated library reference file", required=True)

    parser_fastq_plot_alignments.add_argument("--out_dir", help="Output directory for plots", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_plot_alignments.add_argument("--plot_suf", default=".pdf", help="Plot file suffix (Default: .pdf)")
    parser_fastq_plot_alignments.add_argument("--show", action="store_true", help="Display plots interactively",default=False)
    
    # count_region():
    parser_fastq_count_region.add_argument("--df_ref", help="Annotated reference library file path", required=True)
    parser_fastq_count_region.add_argument("--align_col", help="Align column name in the annotated reference library", required=True)
    parser_fastq_count_region.add_argument("--id_col", help="ID column name in the annotated reference library", required=True)
    parser_fastq_count_region.add_argument("--fastq_dir", help="Directory containing FASTQ files", required=True)
    parser_fastq_count_region.add_argument("--df_motif5", help="5' motif file path", required=True)
    parser_fastq_count_region.add_argument("--df_motif3", help="3' motif file path", required=True)

    parser_fastq_count_region.add_argument("--out_dir", help="Output directory", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_count_region.add_argument("--fastq_col", help="Fastq column name in the annotated reference library (Default: None)", default=None)
    parser_fastq_count_region.add_argument("--match_score", type=float, default=2, help="Score for matches (Default: 2)")
    parser_fastq_count_region.add_argument("--mismatch_score", type=float, default=-1, help="Score for mismatches (Default: -1)")
    parser_fastq_count_region.add_argument("--open_gap_score", type=float, default=-10, help="Gap opening score (Default: -10)")
    parser_fastq_count_region.add_argument("--extend_gap_score", type=float, default=-0.1, help="Gap extension score (Default: -0.1)")
    parser_fastq_count_region.add_argument("--align_dims", type=parse_tuple_int, default=(0, 0), help="Alignment range formatted as 'start,end' (Default: 0,0 = all reads)")
    parser_fastq_count_region.add_argument("--align_ckpt", type=int, default=10000, help="Checkpoint frequency (Default: 10000)")
    parser_fastq_count_region.add_argument("--plot_suf", type=str, help="Plot suffix type (e.g. '.pdf')")
    parser_fastq_count_region.add_argument("--show", action="store_true", help="Display plots interactively", default=False)
    
    # count_alignments():
    parser_fastq_count_alignments.add_argument("--df_ref", help="Annotated reference library file path", required=True)
    parser_fastq_count_alignments.add_argument("--align_col", help="Align column name in the annotated reference library", required=True)
    parser_fastq_count_alignments.add_argument("--id_col", help="ID column name in the annotated reference library", required=True)
    parser_fastq_count_alignments.add_argument("--fastq_dir", help="Directory containing FASTQ files", required=True)

    parser_fastq_count_alignments.add_argument("--out_dir", help="Output directory", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_count_alignments.add_argument("--fastq_col", help="Fastq column name in the annotated reference library (Default: None)", default=None)
    parser_fastq_count_alignments.add_argument("--match_score", type=float, default=2, help="Match score (Default: 2)")
    parser_fastq_count_alignments.add_argument("--mismatch_score", type=float, default=-1, help="Mismatch penalty (Default: -1)")
    parser_fastq_count_alignments.add_argument("--open_gap_score", type=float, default=-10, help="Gap open penalty (Default: -10)")
    parser_fastq_count_alignments.add_argument("--extend_gap_score", type=float, default=-0.1, help="Gap extension penalty (Default: -0.1)")
    parser_fastq_count_alignments.add_argument("--align_dims", type=parse_tuple_int, default=(0, 0), help="Alignment range as 'start,end' (Default: 0,0 = all reads)")
    parser_fastq_count_alignments.add_argument("--align_ckpt", type=int, default=10000, help="Checkpoint frequency for saving alignment progress")
    parser_fastq_count_alignments.add_argument("--plot_suf", type=str, help="Plot file suffix (e.g. .pdf, .png)")
    parser_fastq_count_alignments.add_argument("--show", action="store_true", help="Show plots interactively")
    
    # plot_paired():
    parser_fastq_plot_paired.add_argument("--df", help="Paired region file path", required=True)
    parser_fastq_plot_paired.add_argument("--title", help="Plot title and output filename (without extension)", required=True)
    
    parser_fastq_plot_paired.add_argument("--out_dir", help="Output directory", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_plot_paired.add_argument("--id_col", default="ID", help="Column name for ID (Default: 'ID')")
    parser_fastq_plot_paired.add_argument("--desired_col", default="desired", help="Column name for desired sequences (Default: 'desired')")
    parser_fastq_plot_paired.add_argument("--y", default="count", help="y axis for plots (Default: 'count'; Options: 'count' & 'fraction')")
    parser_fastq_plot_paired.add_argument("--plot_suf", default=".pdf", help="Plot file suffix (e.g. .pdf or .png)")
    parser_fastq_plot_paired.add_argument("--show", action="store_true", help="Display plots interactively")

    # paired_regions():
    parser_fastq_paired_regions.add_argument("--meta_dir", help="Directory containing meta files", required=True)
    parser_fastq_paired_regions.add_argument("--region1_dir", help="Directory with region 1 alignment files", required=True)
    parser_fastq_paired_regions.add_argument("--region2_dir", help="Directory with region 2 alignment files", required=True)

    parser_fastq_paired_regions.add_argument("--out_dir", help="Output directory", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_paired_regions.add_argument("--id_col", default="ID", help="Column name for unique identifiers (Default: 'ID')")
    parser_fastq_paired_regions.add_argument("--desired_col", default="desired", help="Column name for desired sequences (Default: 'desired')")
    parser_fastq_paired_regions.add_argument("--region1_alignment_col", default="r1_alignment", help="Column name for region 1 alignment data")
    parser_fastq_paired_regions.add_argument("--region2_alignment_col", default="r2_alignment", help="Column name for region 2 alignment data")
    parser_fastq_paired_regions.add_argument("--reads_aligned_col", default="reads_aligned", help="Column name for aligned reads (Default: 'reads_aligned')")
    parser_fastq_paired_regions.add_argument("--reads_processed_col", default="reads_processed", help="Column name for processed reads (Default: 'reads_processed')")
    parser_fastq_paired_regions.add_argument("--plot_suf", default=".pdf", help="Plot file suffix (e.g., .pdf, .png)")
    parser_fastq_paired_regions.add_argument("--show", action="store_true", help="Display plots interactively")
    parser_fastq_paired_regions.add_argument("--return_dc", action="store_true", help="Return paired/unpaired dataframe")

    # count_signatures():
    parser_fastq_count_signatures.add_argument("--df_ref", help="Annotated reference library file path", required=True)
    parser_fastq_count_signatures.add_argument("--fastq_dir", help="Directory containing FASTQ files", required=True)
    
    parser_fastq_count_signatures_group = parser_fastq_count_signatures.add_mutually_exclusive_group(required=True)
    parser_fastq_count_signatures_group.add_argument("--target_sequence", help="[Required (Option 1)] Target sequence; retrieved from input file if not provided")
    parser_fastq_count_signatures_group.add_argument("--in_file", help="[Required (Option 2)] Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence,aa_index (column names required)")
    
    parser_fastq_count_signatures.add_argument("--df_motif5", help="5' motif file path", default=argparse.SUPPRESS)
    parser_fastq_count_signatures.add_argument("--df_motif3", help="3' motif file path", default=argparse.SUPPRESS)
    parser_fastq_count_signatures.add_argument("--meta", help="Meta file path", default=argparse.SUPPRESS)
    parser_fastq_count_signatures.add_argument("--out_dir", help="Output directory", default='../out/')
    parser_fastq_count_signatures.add_argument("--out_file", help="Output filename", default=f"{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_count_signatures.csv")
    parser_fastq_count_signatures.add_argument("--signature_col", help="Signature column name in the annotated reference library (Default: 'Signature')", default='Signature')
    parser_fastq_count_signatures.add_argument("--id_col", help="ID column name in the annotated reference library (Default: 'ID')", default='ID')
    parser_fastq_count_signatures.add_argument("--edit_col", help="Edit column name in the annotated reference library (Default: 'Edit')", default='Edit')
    parser_fastq_count_signatures.add_argument("--fastq_col", help="Fastq column name in the annotated reference library (Default: None)", default=None)
    parser_fastq_count_signatures.add_argument("--match_score", type=float, default=2, help="Score for matches (Default: 2)")
    parser_fastq_count_signatures.add_argument("--mismatch_score", type=float, default=-1, help="Score for mismatches (Default: -1)")
    parser_fastq_count_signatures.add_argument("--open_gap_score", type=float, default=-10, help="Gap opening score (Default: -10)")
    parser_fastq_count_signatures.add_argument("--extend_gap_score", type=float, default=-0.1, help="Gap extension score (Default: -0.1)")
    parser_fastq_count_signatures.add_argument("--align_dims", type=parse_tuple_int, default=(0, 0), help="Alignment range formatted as 'start,end' (Default: 0,0 = all reads)")
    parser_fastq_count_signatures.add_argument("--align_ckpt", type=int, default=10000, help="Checkpoint frequency (Default: 10000)")
    parser_fastq_count_signatures.add_argument("--save_alignments", action="store_true", help="Save alignments (Default: False, save memory)", default=False)
    parser_fastq_count_signatures.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)
    parser_fastq_count_signatures.add_argument("--plot_suf", type=str, help="Plot suffix type (Default: 'pdf')", default='.pdf')
    parser_fastq_count_signatures.add_argument("--show", action="store_true", help="Display plots interactively", default=False)
    
    # editing_per_library():
    parser_fastq_editing_per_library.add_argument("--edit_dc", help="Path to directory with edit outcomes files", required=True)
    parser_fastq_editing_per_library.add_argument("--paired_regions_dc", help="Path to directory with paired regions files", required=True)
    parser_fastq_editing_per_library.add_argument("--fastq_ids", help="Path to file containing fastq IDs for 'genotyping' and 'paired_regions'", required=True)

    parser_fastq_editing_per_library.add_argument("--out_dir", type=str, help="Output directory to save results (Default: ../out/date_time)", default=f"../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}")
    parser_fastq_editing_per_library.add_argument("--count", default="count", help="Column to use for epeg-ngRNA counts (Default: 'count')")
    parser_fastq_editing_per_library.add_argument("--psuedocount", type=int, default=1, help="Pseudocount to add to all counts (Default: 1)")

    # extract_umis():
    parser_fastq_extract_umis.add_argument("--fastq_dir", help="Directory containing FASTQ files", required=True)

    parser_fastq_extract_umis.add_argument("--out_dir", help="Output directory (Default: ./extract_umis)", default=f'./extract_umis')
    parser_fastq_extract_umis.add_argument("--bc_pattern", help="UMI barcode pattern (Default: NNNNNNNNNNNNNNNN)", default="NNNNNNNNNNNNNNNN")
    parser_fastq_extract_umis.add_argument("--env", help="Conda environment with umi_tools installed (Default: umi_tools)", default="umi_tools")

    # trim_motifs():
    parser_fastq_trim_motifs.add_argument("--fastq_dir", help="Directory containing FASTQ files (with UMIs extracted)", required=True)

    parser_fastq_trim_motifs.add_argument("--out_dir", help="Output directory (Default: ./trim_motifs)", default=f'./trim_motifs')

    parser_fastq_trim_motifs.add_argument("--in_file", help="[Required option 1] Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)", default=argparse.SUPPRESS)
    parser_fastq_trim_motifs.add_argument("--motif5", help="[Required option 2] 5' motif sequence to trim", default=argparse.SUPPRESS)
    parser_fastq_trim_motifs.add_argument("--motif3", help="[Required option 2] 3' motif sequence to trim", default=argparse.SUPPRESS)

    parser_fastq_trim_motifs.add_argument("--motif_length", type=int, help="Trim 'in_file' motifs to this length (Default: 21)", default=21)
    parser_fastq_trim_motifs.add_argument("--error_rate", type=float, help="Maximum error rate allowed in each motif (Default: 0.1 = 10%%)", default=0.1)
    parser_fastq_trim_motifs.add_argument("--max_distance", type=int, help="Maximum number of errors allowed in each motif (Default: 2)", default=2)
    parser_fastq_trim_motifs.add_argument("--env", help="Conda environment with cutadapt installed (Default: cutadapt)", default="umi_tools")

    # make_sams():
    parser_fastq_make_sams.add_argument("--fastq_dir", help="Directory containing FASTQ files", required=True)

    parser_fastq_make_sams.add_argument("--out_dir", help="Output directory (Default: ./make_sams)", default=f'./make_sams')

    parser_fastq_make_sams_group = parser_fastq_make_sams.add_mutually_exclusive_group(required=True)
    parser_fastq_make_sams_group.add_argument("--in_file", help="Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)", default=argparse.SUPPRESS)
    parser_fastq_make_sams_group.add_argument("--fasta", help="Reference FASTA file for alignment", default=argparse.SUPPRESS)

    parser_fastq_make_sams.add_argument("--sensitivity", choices=["very-sensitive", "sensitive", "fast", "very-fast", "very-sensitive-local", "sensitive-local", "fast-local", "very-fast-local"], 
                                        default="very-sensitive", help="Bowtie2 sensitivity setting (Default: very-sensitive)")

    parser_fastq_make_sams.add_argument("--env", help="Conda environment with bowtie2 installed (Default: umi_tools)", default="umi_tools")

    # make_bams():
    parser_fastq_make_bams.add_argument("--sam_dir", help="Directory containing SAM files", required=True)

    parser_fastq_make_bams.add_argument("--out_dir", help="Output directory (Default: ./make_bams)", default=f'./make_bams')
    parser_fastq_make_bams.add_argument("--env", help="Conda environment with samtools installed (Default: umi_tools)", default="umi_tools")

    # bam_umi_tags():
    parser_fastq_bam_umi_tags.add_argument("--bam_dir", help="Directory containing bam files", required=True)

    parser_fastq_bam_umi_tags.add_argument("--out_dir", help="Output directory (Default: ./bam_umi_tags)", default=f'./bam_umi_tags')
    parser_fastq_bam_umi_tags.add_argument("--env", help="Conda environment with fgbio installed (Default: umi_tools)", default="umi_tools")

    # group_umis():
    parser_fastq_group_umis.add_argument("--bam_dir", help="Directory containing BAM files", required=True)

    parser_fastq_group_umis.add_argument("--out_dir", help="Output directory (Default: ./group_umis)", default=f'./group_umis')
    parser_fastq_group_umis.add_argument("--strategy", choices=["Identical","Edit","Adjacency", "Paired"], help="umi grouping strategy (Default: Adjacency)", default="Adjacency")
    parser_fastq_group_umis.add_argument("--edits", type=int, help="Maximum edit distance to group UMIs (Default: 1)", default=1)
    parser_fastq_group_umis.add_argument("--env", help="Conda environment with fgbio installed (Default: umi_tools)", default="umi_tools")

    # consensus_umis():
    parser_fastq_consensus_umis.add_argument("--bam_dir", help="Directory containing grouped BAM files", required=True)

    parser_fastq_consensus_umis.add_argument("--out_dir", help="Output directory (Default: ./consensus_umis)", default=f'./consensus_umis')
    parser_fastq_consensus_umis.add_argument("--min_reads", type=int, help="Minimum reads per UMI to call consensus (Default: 1)", default=1)
    parser_fastq_consensus_umis.add_argument("--env", help="Conda environment with fgbio installed (Default: umi_tools)", default="umi_tools")

    # bam_to_fastq():
    parser_fastq_bam_to_fastq.add_argument("--bam_dir", help="Directory containing BAM files", required=True)
    parser_fastq_bam_to_fastq.add_argument("--out_dir", help="Output directory (Default: ./bam_to_fastq)", default=f'./bam_to_fastq')
    parser_fastq_bam_to_fastq.add_argument("--env", help="Conda environment with samtools installed (Default: umi_tools)", default="umi_tools")
    
    # Set defaults
    parser_fastq_revcom.set_defaults(func=fq.revcom_fastqs)
    parser_fastq_unzip.set_defaults(func=fq.unzip_fastqs)
    parser_fastq_comb.set_defaults(func=fq.comb_fastqs)
    parser_fastq_genotyping.set_defaults(func=fq.genotyping)
    parser_fastq_abundances.set_defaults(func=fq.abundances)
    parser_fastq_count_motif.set_defaults(func=fq.count_motif)
    parser_fastq_plot_motif.set_defaults(func=fq.plot_motif)
    parser_fastq_plot_alignments.set_defaults(func=fq.plot_alignments)
    parser_fastq_count_region.set_defaults(func=fq.count_region)
    parser_fastq_count_alignments.set_defaults(func=fq.count_alignments)
    parser_fastq_plot_paired.set_defaults(func=fq.plot_paired)
    parser_fastq_paired_regions.set_defaults(func=fq.paired_regions)
    parser_fastq_count_signatures.set_defaults(func=fq.count_signatures)
    parser_fastq_editing_per_library.set_defaults(func=fq.editing_per_library)
    parser_fastq_extract_umis.set_defaults(func=fq.extract_umis)
    parser_fastq_trim_motifs.set_defaults(func=fq.trim_motifs)
    parser_fastq_make_sams.set_defaults(func=fq.make_sams)
    parser_fastq_make_bams.set_defaults(func=fq.make_bams)
    parser_fastq_bam_umi_tags.set_defaults(func=fq.bam_umi_tags)
    parser_fastq_group_umis.set_defaults(func=fq.group_umis)
    parser_fastq_consensus_umis.set_defaults(func=fq.consensus_umis)
    parser_fastq_bam_to_fastq.set_defaults(func=fq.bam_to_fastq)

    '''
    Add pe.py
    - prime_designer(): Execute PrimeDesign saturation mutagenesis for EDMS using Docker (NEED TO BE RUNNING DESKTOP APP)
    - pilot_screen(): Create pilot screen for EDMS
    - epegRNA_linkers(): Generate epegRNA linkers between PBS and 3' hairpin motif & finish annotations
    - merge(): rejoins epeg/ngRNAs & creates ngRNA_groups
    - sensor_designer(): design pegRNA sensors
    - rtt_designer(): design all possible RTT for given spacer & PBS (WT, single insertions, & single deletions)
    - pegRNA_outcome(): confirm that pegRNAs should create the predicted edit
    - pegRNA_signature(): create signatures for pegRNA outcomes using alignments
    '''
    parser_pe = subparsers.add_parser("pe", help="Prime Editing")
    subparsers_pe = parser_pe.add_subparsers()

    parser_pe_prime_designer = subparsers_pe.add_parser("prime_designer", help="Execute PrimeDesign saturation mutagenesis for EDMS")
    parser_pe_pilot_screen = subparsers_pe.add_parser("pilot_screen", help="Determine pilot screen for EDMS")
    parser_pe_epegRNA_linkers = subparsers_pe.add_parser("epegRNA_linkers", help="Generate epegRNA linkers between PBS and 3' hairpin motif")
    parser_pe_merge = subparsers_pe.add_parser("merge", help="rejoins epeg/ngRNAs & creates ngRNA groups")
    parser_pe_sensor_designer = subparsers_pe.add_parser("sensor_designer", help='Design pegRNA sensors')
    parser_pe_rtt_designer = subparsers_pe.add_parser("rtt_designer", help="Design all possible RTT for given spacer & PBS (WT, single insertions, & single deletions)")
    parser_pe_pegRNA_outcome = subparsers_pe.add_parser("pegRNA_outcome", help="Confirm that pegRNAs should create the predicted edit")
    parser_pe_pegRNA_signature = subparsers_pe.add_parser("pegRNA_signature", help="Create signatures for pegRNA outcomes using alignments")

    # prime_designer():
    parser_pe_prime_designer.add_argument("--name", type=str, dest='target_name',help="Name of the target", required=True)
    parser_pe_prime_designer.add_argument("--flank5", type=str, dest='flank5_sequence', help="5' flank sequence (in-frame, length divisible by 3)", required=True)
    parser_pe_prime_designer.add_argument("--target", type=str, dest='target_sequence', help="Target sequence (in-frame, length divisible by 3)", required=True)
    parser_pe_prime_designer.add_argument("--flank3", type=str, dest='flank3_sequence', help="3' flank sequence (in-frame, length divisible by 3)", required=True)

    parser_pe_prime_designer.add_argument("--pbs_lengths", type=int, dest='pbs_length_pooled_ls', nargs="+", default=[11,13,15],
                        help="List of PBS lengths (Default: 11,13,15)")
    parser_pe_prime_designer.add_argument("--rtt_max_length", type=int, dest='rtt_max_length_pooled', default=50,
                        help="Maximum RTT length to design pegRNAs (Default: 50 nt)")
    parser_pe_prime_designer.add_argument("--no_silent_mutation", dest="silent_mutation", action="store_false",
                        help="Disable silent mutation")
    parser_pe_prime_designer.add_argument("--number_of_pegrnas", type=int, default=1,
                        help="Max number of pegRNAs to design (Default: 1)")
    parser_pe_prime_designer.add_argument("--number_of_ngrnas", type=int, default=3,
                        help="Max number of ngRNAs to design (Default: 3)")
    parser_pe_prime_designer.add_argument("--scaffold_sequence", type=str, default="GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
                        help="sgRNA scaffold sequence (Default: SpCas9 flip + extend")
    parser_pe_prime_designer.add_argument("--aa_index", type=int, default=1,
                        help="Index of 1st amino acid in target sequence (Default: 1)")
    parser_pe_prime_designer.add_argument("--enzymes", type=str, nargs="+", help="list of type IIS RE enzymes (i.e., Esp3I, BsaI, BspMI) to check for in pegRNAs and ngRNAs (Default: ['Esp3I'])", default=['Esp3I'])
    parser_pe_prime_designer.add_argument("--dont_replace", action='store_false',dest='replace', help="Do not replace pegRNAs and remove ngRNAs with RE enzyme sites", default=True)

    # Pilot_Screen():
    parser_pe_pilot_screen.add_argument("--pegRNAs", type=str, dest='pegRNAs_dir',help="Directory with pegRNAs from prime_designer() output", required=True)
    parser_pe_pilot_screen.add_argument("--mutations", type=str, dest='mutations_pt', help="Path to mutations file (COSMIC or ClinVar)", required=True)
    
    parser_pe_pilot_screen.add_argument("--database", type=str, choices=['COSMIC', 'ClinVar'], default='COSMIC', help="Database to use for priority mutations (Default: 'COSMIC')")
    parser_pe_pilot_screen.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    # epegRNA_linkers():
    parser_pe_epegRNA_linkers.add_argument('--pegRNAs', help='Path to pegRNAs file',required=True)

    parser_pe_epegRNA_linkers.add_argument('--epegRNA_motif_sequence', default='CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA', help='epegRNA motif sequence (Default: tevopreQ1)')
    parser_pe_epegRNA_linkers.add_argument('--linker_pattern', type=str, default=argparse.SUPPRESS, help='epegRNA linker pattern (Default: NNNNNNNN)')
    parser_pe_epegRNA_linkers.add_argument('--excluded_motifs', type=str, nargs="+", default=['Esp3I'], help="list of motifs or type IIS RE enzymes (i.e., Esp3I, BsaI, BspMI) to exclude from linker generation (Default: ['Esp3I'])")
    parser_pe_epegRNA_linkers.add_argument('--ckpt_dir', type=str, help='Checkpoint directory path (Default: ../epegRNAs/ckpt)', default='../epegRNAs/ckpt')
    parser_pe_epegRNA_linkers.add_argument('--ckpt_file', help='Checkpoint file name (Default: YYMMDD_HHMMSS_epegRNA_linkers.csv)', default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_epegRNA_linkers.csv')
    parser_pe_epegRNA_linkers.add_argument('--ckpt_pt', type=str, default='', help='Previous checkpoint full path (Example: ../epegRNAs/ckpt/YYMMDD_HHMMSS_epegRNA_linkers.csv)')
    parser_pe_epegRNA_linkers.add_argument("--out_dir", type=str, help="Output directory (Default: ../epegRNAs)", default='../epegRNAs')
    parser_pe_epegRNA_linkers.add_argument("--out_file", type=str, help="Name of the output file (Default: epegRNAs.csv)", default='epegRNAs.csv')
    parser_pe_epegRNA_linkers.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    # merge():
    parser_pe_merge.add_argument("--epegRNAs", type=str, help="Directory or file with epegRNAs", required=True)
    parser_pe_merge.add_argument("--ngRNAs", type=str, help="Directory or file with ngRNAs", required=True)
    
    parser_pe_merge.add_argument("--ngRNAs_groups", type=str, dest='ngRNAs_groups_max', help="Maximum # of ngRNAs per epegRNA (Default: 3)", default=3)
    parser_pe_merge.add_argument("--epegRNA_suffix", type=str, help="Suffix for epegRNAs columns (Default: _epegRNA)", default='_epegRNA')
    parser_pe_merge.add_argument("--ngRNA_suffix", type=str, help="Suffix for ngRNAs columns (Default: _ngRNA)", default='_ngRNA')
    parser_pe_merge.add_argument("--out_dir", type=str, dest='dir', help="Output directory (Default: ../epeg_ngRNAs)", default='../epeg_ngRNAs')
    parser_pe_merge.add_argument("--out_file", type=str, dest='file', help="Name of the output file (Default: epeg_ngRNAs.csv)", default='epeg_ngRNAs.csv')
    parser_pe_merge.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    # sensor_designer():
    parser_pe_sensor_designer.add_argument("--pegRNAs", type=str, help="Path to pegRNAs file", required=True)
    parser_pe_sensor_designer.add_argument("--in_file", type=str, help="Path to PrimeDesign input file (required columns: target_name, target_sequence, aa_index)", required=True)

    parser_pe_sensor_designer.add_argument("--sensor_length", type=int, default=60, help="Total length of the sensor in bp (Default: 60)")
    parser_pe_sensor_designer.add_argument("--before_spacer", type=int, default=5, help="Amount of nucleotide context to put before the protospacer in the sensor (Default = 5)")
    parser_pe_sensor_designer.add_argument("--sensor_orientation", type=str, default='revcom', help="Orientation of the sensor relative to the protospacer (Options: 'revcom' [Default b/c minimize recombination] or ’forward’")
    parser_pe_sensor_designer.add_argument("--out_dir", type=str, help="Output directory (Default: ../pegRNAs_tester)", default='../sensor_designer')
    parser_pe_sensor_designer.add_argument("--out_file", type=str, help="Name of the output file (Default: pegRNAs.csv)", default='pegRNAs.csv')
    parser_pe_sensor_designer.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    # rtt_designer():
    parser_pe_rtt_designer.add_argument("--pegRNAs", type=str, help="Path to pegRNAs file", required=True)
    parser_pe_rtt_designer.add_argument("--in_file", type=str, help="Path to PrimeDesign input file (required columns: target_name, target_sequence, aa_index)", required=True)
    
    parser_pe_rtt_designer.add_argument("--RTT_length", type=int, default=39, help="Length of RTT (Default: 39)")
    parser_pe_rtt_designer.add_argument("--include_WT", action='store_true', default=False, help="include wildtype RTT (Default: False)")
    parser_pe_rtt_designer.add_argument("--out_dir", type=str, help="Output directory (Default: ../RTT_Designer)", default='../RTT_Designer')
    parser_pe_rtt_designer.add_argument("--out_file", type=str, help="Name of the output file (Default: pegRNAs.csv)", default='pegRNAs.csv')
    parser_pe_rtt_designer.add_argument("--comments", action='store_true', help="Print comments (Default: False)", default=False)
    parser_pe_rtt_designer.add_argument("--enzymes", type=str, nargs="+", help="list of type IIS RE enzymes (i.e., Esp3I, BsaI, BspMI) to check for in pegRNAs and ngRNAs (Default: ['Esp3I'])", default=['Esp3I'])
    parser_pe_rtt_designer.add_argument("--dont_replace", action='store_false',dest='replace', help="Do not replace pegRNAs with RE enzyme sites", default=True)
    parser_pe_rtt_designer.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    # pegRNA_outcome():
    parser_pe_pegRNA_outcome.add_argument("--pegRNAs", type=str, help="Path to pegRNAs file", required=True)
    parser_pe_pegRNA_outcome.add_argument("--in_file", type=str, help="Path to PrimeDesign input file (required columns: target_name, target_sequence, aa_index)", required=True)

    parser_pe_pegRNA_outcome.add_argument("--out_dir", type=str, help="Output directory (Default: ../pegRNA_outcome)", default='../pegRNA_outcome')
    parser_pe_pegRNA_outcome.add_argument("--out_file", type=str, help="Name of the output file (Default: pegRNAs.csv)", default='pegRNAs.csv')
    parser_pe_pegRNA_outcome.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    parser_pe_pegRNA_outcome.add_argument("--match_score", type=float, help="Match score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_outcome.add_argument("--mismatch_score", type=float, help="Mismatch score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_outcome.add_argument("--open_gap_score", type=float, help="Open gap score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_outcome.add_argument("--extend_gap_score", type=float, help="Extend gap score for pairwise alignment", default=argparse.SUPPRESS)

    # pegRNA_signature():
    parser_pe_pegRNA_signature.add_argument("--pegRNAs", type=str, help="Path to pegRNAs file", required=True)
    parser_pe_pegRNA_signature.add_argument("--in_file", type=str, help="Path to PrimeDesign input file (required columns: target_name, target_sequence, aa_index)", required=True)

    parser_pe_pegRNA_signature.add_argument("--out_dir", type=str, help="Output directory (Default: ../pegRNA_signature)", default='../pegRNA_signature')
    parser_pe_pegRNA_signature.add_argument("--out_file", type=str, help="Name of the output file (Default: pegRNAs.csv)", default='pegRNAs.csv')
    parser_pe_pegRNA_signature.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    parser_pe_pegRNA_signature.add_argument("--match_score", type=float, help="Match score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_signature.add_argument("--mismatch_score", type=float, help="Mismatch score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_signature.add_argument("--open_gap_score", type=float, help="Open gap score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_signature.add_argument("--extend_gap_score", type=float, help="Extend gap score for pairwise alignment", default=argparse.SUPPRESS)

    # Set defaults
    parser_pe_prime_designer.set_defaults(func=pe.prime_designer)
    parser_pe_pilot_screen.set_defaults(func=pe.pilot_screen)
    parser_pe_epegRNA_linkers.set_defaults(func=pe.epegRNA_linkers)
    parser_pe_merge.set_defaults(func=pe.merge)
    parser_pe_sensor_designer.set_defaults(func=pe.sensor_designer)
    parser_pe_rtt_designer.set_defaults(func=pe.rtt_designer)
    parser_pe_pegRNA_outcome.set_defaults(func=pe.pegRNA_outcome)
    parser_pe_pegRNA_signature.set_defaults(func=pe.pegRNA_signature)

    # Enable autocomplete
    argcomplete.autocomplete(parser)

    # Parse all arguments
    args = parser.parse_args()
    if args.command == 'autocomplete': # Run autocomplete script
        sys.exit(utils.run_bundled_script())
    args_dict = vars(args)
    func = args_dict.pop("func")
    args_dict.pop("command", None)  # Remove 'command' if it exists (required for autocomplete)
    func(**args_dict)    