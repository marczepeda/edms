'''
Module: main.py
Author: Marc Zepeda
Created: 2025-04-12
Description: Endogenous Deep Mutational Scans

Usage:
[Supporting methods]
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
import argparse
import ast
import datetime

from . import config

from .gen import plot as p
from .gen import stat as st
from .gen import io
from .gen import cli

from .dat import cosmic as co 
from .dat import cvar

from .bio import ngs
from .bio import transfect as tf
from .bio import qPCR
from .bio import clone as cl
from .bio import fastq as fq
from .bio import pe

# Supporting methods
'''
    add_common_plot_scat_args(subparser): Add common arguments for scatter plot related graphs
'''
def add_common_plot_scat_args(subparser):

    # scat(): Required arguments
    subparser.add_argument("--df", help="Input file", type=str, required=True)
    subparser.add_argument("--x", help="X-axis column", type=str, required=True)
    subparser.add_argument("--y", help="Y-axis column", type=str, required=True)

    # Optional core arguments
    subparser.add_argument("--cols", type=str, help="Color column name")
    subparser.add_argument("--cols_ord", nargs="+", help="Column order (list)")
    subparser.add_argument("--cols_exclude", nargs="+", help="Columns to exclude")
    subparser.add_argument("--stys", type=str, help="Style column name")

    subparser.add_argument("--dir", help="Output directory path", type=str, default='.')
    subparser.add_argument("--file", help="Output file name", type=str, required=False, default='plot_scat.png')
    subparser.add_argument("--palette_or_cmap", type=str, default="colorblind", help="Color palette or colormap")
    subparser.add_argument("--edgecol", type=str, default="black", help="Edge color")

    # Figure appearance
    subparser.add_argument("--figsize", nargs=2, type=int, default=(10,6), help="Figure size (width height)")
    subparser.add_argument("--title", type=str, default="", help="Plot title")
    subparser.add_argument("--title_size", type=int, default=18)
    subparser.add_argument("--title_weight", type=str, default="bold")
    subparser.add_argument("--title_font", type=str, default="Arial")

    # X-axis settings
    subparser.add_argument("--x_axis", type=str, default="", help="X-axis label")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_axis_font", type=str, default="Arial")
    subparser.add_argument("--x_axis_scale", type=str, default="linear")
    subparser.add_argument("--x_axis_dims", nargs=2, type=float, default=(0,0))
    subparser.add_argument("--x_ticks_rot", type=int, default=0)
    subparser.add_argument("--x_ticks_font", type=str, default="Arial")
    subparser.add_argument("--x_ticks", nargs="+", help="X-axis ticks")

    # Y-axis settings
    subparser.add_argument("--y_axis", type=str, default="", help="Y-axis label")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_axis_font", type=str, default="Arial")
    subparser.add_argument("--y_axis_scale", type=str, default="linear")
    subparser.add_argument("--y_axis_dims", nargs=2, type=float, default=(0,0))
    subparser.add_argument("--y_ticks_rot", type=int, default=0)
    subparser.add_argument("--y_ticks_font", type=str, default="Arial")
    subparser.add_argument("--y_ticks", nargs="+", help="Y-axis ticks")

    # Legend settings
    subparser.add_argument("--legend_title", type=str, default="")
    subparser.add_argument("--legend_title_size", type=int, default=12)
    subparser.add_argument("--legend_size", type=int, default=9)
    subparser.add_argument("--legend_bbox_to_anchor", nargs=2, type=float, default=(1,1))
    subparser.add_argument("--legend_loc", type=str, default="upper left")
    subparser.add_argument("--legend_items", nargs=2, type=int, default=(0,0))
    subparser.add_argument("--legend_ncol", type=int, default=1)

    # Display and formatting
    subparser.add_argument("--show", action="store_true", help="Show the plot")
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize legend/label names with spaces")

'''
    add_common_plot_cat_args(subparser): Add common arguments for category dependent graphs
'''
def add_common_plot_cat_args(subparser):
    
    # cat(): Required arguments
    subparser.add_argument("--df", help="Input file", type=str, required=True)

    # Optional core arguments
    subparser.add_argument("--x", help="X-axis column", type=str, default="")
    subparser.add_argument("--y", help="Y-axis column", type=str, default="")
    subparser.add_argument("--cols", type=str, help="Column used for color grouping")
    subparser.add_argument("--cols_ord", nargs="+", help="Custom order for color values")
    subparser.add_argument("--cols_exclude", nargs="+", help="Values to exclude from color column")

    subparser.add_argument("--file", type=str, help="Output filename")
    subparser.add_argument("--dir", type=str, help="Output directory path")
    subparser.add_argument("--palette_or_cmap", type=str, default="colorblind", help="Color palette or colormap")
    subparser.add_argument("--edgecol", type=str, default="black", help="Edge color")

    # Error bar and style options
    subparser.add_argument("--lw", type=int, default=1, help="Line width")
    subparser.add_argument("--errorbar", type=str, default="sd", help="Error bar type (e.g., sd)")
    subparser.add_argument("--errwid", type=float, default=1, help="Error bar width")
    subparser.add_argument("--errcap", type=float, default=0.1, help="Error bar cap size")

    # Figure appearance
    subparser.add_argument("--figsize", nargs=2, type=int, default=(10, 6), help="Figure size (width height)")
    subparser.add_argument("--title", type=str, default="")
    subparser.add_argument("--title_size", type=int, default=18)
    subparser.add_argument("--title_weight", type=str, default="bold")
    subparser.add_argument("--title_font", type=str, default="Arial")

    # X-axis settings
    subparser.add_argument("--x_axis", type=str, default="")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_axis_font", type=str, default="Arial")
    subparser.add_argument("--x_axis_scale", type=str, default="linear")
    subparser.add_argument("--x_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--x_ticks_rot", type=int, default=0)
    subparser.add_argument("--x_ticks_font", type=str, default="Arial")
    subparser.add_argument("--x_ticks", nargs="+", help="X-axis ticks")

    # Y-axis settings
    subparser.add_argument("--y_axis", type=str, default="")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_axis_font", type=str, default="Arial")
    subparser.add_argument("--y_axis_scale", type=str, default="linear")
    subparser.add_argument("--y_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--y_ticks_rot", type=int, default=0)
    subparser.add_argument("--y_ticks_font", type=str, default="Arial")
    subparser.add_argument("--y_ticks", nargs="+", help="Y-axis ticks")

    # Legend settings
    subparser.add_argument("--legend_title", type=str, default="")
    subparser.add_argument("--legend_title_size", type=int, default=12)
    subparser.add_argument("--legend_size", type=int, default=9)
    subparser.add_argument("--legend_bbox_to_anchor", nargs=2, type=float, default=(1, 1))
    subparser.add_argument("--legend_loc", type=str, default="upper left")
    subparser.add_argument("--legend_items", nargs=2, type=int, default=(0, 0))
    subparser.add_argument("--legend_ncol", type=int, default=1)

    # Display and formatting
    subparser.add_argument("--show", action="store_true", help="Show the plot")
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize labels/legends")

'''
    add_common_plot_dist_args(subparser): Add common arguments for distribution graphs
'''
def add_common_plot_dist_args(subparser):
    ''''''
    # dist(): Required argument
    subparser.add_argument("--df", help="Input file", type=str, required=True)
    subparser.add_argument("--x", type=str, help="X-axis column", required=True)

    # Optional core arguments
    subparser.add_argument("--cols", type=str, help="Color column")
    subparser.add_argument("--cols_ord", nargs="+", help="Custom color order")
    subparser.add_argument("--cols_exclude", nargs="+", help="Values to exclude from color column")

    # Plot customization
    subparser.add_argument("--bins", type=int, default=40, help="Number of bins for histogram")
    subparser.add_argument("--log10_low", type=int, default=0, help="Log10 lower limit for scale")
    subparser.add_argument("--palette_or_cmap", type=str, default="colorblind")
    subparser.add_argument("--edgecol", type=str, default="black")
    subparser.add_argument("--lw", type=int, default=1, help="Line width")
    subparser.add_argument("--ht", type=float, default=1.5, help="Plot height")
    subparser.add_argument("--asp", type=int, default=5, help="Aspect ratio")
    subparser.add_argument("--tp", type=float, default=0.8, help="Top padding")
    subparser.add_argument("--hs", type=int, default=0, help="Hspace between plots")
    subparser.add_argument("--des", action="store_true", help="Remove plot spines (despine)")

    # File output
    subparser.add_argument("--file", type=str, help="Output file name")
    subparser.add_argument("--dir", type=str, help="Output directory")

    # Figure appearance
    subparser.add_argument("--figsize", nargs=2, type=int, default=(10, 6), help="Figure size (width height)")
    subparser.add_argument("--title", type=str, default="")
    subparser.add_argument("--title_size", type=int, default=18)
    subparser.add_argument("--title_weight", type=str, default="bold")
    subparser.add_argument("--title_font", type=str, default="Arial")

    # X-axis
    subparser.add_argument("--x_axis", type=str, default="")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_axis_font", type=str, default="Arial")
    subparser.add_argument("--x_axis_scale", type=str, default="linear")
    subparser.add_argument("--x_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--x_ticks_rot", type=int, default=0)
    subparser.add_argument("--x_ticks_font", type=str, default="Arial")
    subparser.add_argument("--x_ticks", nargs="+", help="X-axis tick values")

    # Y-axis
    subparser.add_argument("--y_axis", type=str, default="")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_axis_font", type=str, default="Arial")
    subparser.add_argument("--y_axis_scale", type=str, default="linear")
    subparser.add_argument("--y_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--y_ticks_rot", type=int, default=0)
    subparser.add_argument("--y_ticks_font", type=str, default="Arial")
    subparser.add_argument("--y_ticks", nargs="+", help="Y-axis tick values")

    # Legend
    subparser.add_argument("--legend_title", type=str, default="")
    subparser.add_argument("--legend_title_size", type=int, default=12)
    subparser.add_argument("--legend_size", type=int, default=9)
    subparser.add_argument("--legend_bbox_to_anchor", nargs=2, type=float, default=(1, 1))
    subparser.add_argument("--legend_loc", type=str, default="upper left")
    subparser.add_argument("--legend_items", nargs=2, type=int, default=(0, 0))
    subparser.add_argument("--legend_ncol", type=int, default=1)

    # Final display
    subparser.add_argument("--show", action="store_true", help="Show the plot")
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize labels/legend with spacing")

'''
    add_common_plot_heat_args(subparser): Add common arguments for heatmap graphs
'''
def add_common_plot_heat_args(subparser):
    
    # Required arguments
    subparser.add_argument("--df", help="Input file", type=str, required=True)
    
    # Optional arguments
    subparser.add_argument("--x", type=str, help="X-axis column name to split tidy-formatted dataframe into a dictionary pivot-formatted dataframes")
    subparser.add_argument("--y", type=str, help="Y-axis column name to split tidy-formatted dataframe into a dictionary pivot-formatted dataframes")
    subparser.add_argument("--vars", type=str, help="Variable column name to split tidy-formatted dataframe into a dictionary pivot-formatted dataframes")
    subparser.add_argument("--vals", type=str, help="Value column name to split tidy-formatted dataframe into a dictionary pivot-formatted dataframes")
    subparser.add_argument("--vals_dims", nargs=2, type=float, help="Min/max for values (vmin vmax)")

    subparser.add_argument("--file", type=str, help="Output filename")
    subparser.add_argument("--dir", type=str, help="Output directory path")
    subparser.add_argument("--edgecol", type=str, default="black", help="Edge color")
    subparser.add_argument("--lw", type=int, default=1, help="Line width")

    subparser.add_argument("--annot", action="store_true", help="Show value annotations in cells")
    subparser.add_argument("--cmap", type=str, default="Reds", help="Colormap name")
    subparser.add_argument("--sq", action="store_true", help="Use square cells (aspect ratio 1:1)")
    subparser.add_argument("--cbar", action="store_true", help="Show colorbar")

    # Title and size
    subparser.add_argument("--title", type=str, default="")
    subparser.add_argument("--title_size", type=int, default=18)
    subparser.add_argument("--title_weight", type=str, default="bold")
    subparser.add_argument("--title_font", type=str, default="Arial")
    subparser.add_argument("--figsize", nargs=2, type=int, default=(10, 6), help="Figure size (width height)")

    # X-axis
    subparser.add_argument("--x_axis", type=str, default="")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_axis_font", type=str, default="Arial")
    subparser.add_argument("--x_ticks_rot", type=int, default=0)
    subparser.add_argument("--x_ticks_font", type=str, default="Arial")

    # Y-axis
    subparser.add_argument("--y_axis", type=str, default="")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_axis_font", type=str, default="Arial")
    subparser.add_argument("--y_ticks_rot", type=int, default=0)
    subparser.add_argument("--y_ticks_font", type=str, default="Arial")

    # Final display
    subparser.add_argument("--show", action="store_true", help="Show the plot")
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize labels/legend with spacing")

'''
    add_common_plot_stack_args(subparser): Add common arguments for stacked bar plot
'''
def add_common_plot_stack_args(subparser):

    # Required arguments
    subparser.add_argument("--df", type=str, help="Input file", required=True)
    subparser.add_argument("--x", type=str, help="X-axis column")
    subparser.add_argument("--y", type=str, help="Y-axis column")
    subparser.add_argument("--cols", type=str, help="Color column")

    # Optional parameters
    subparser.add_argument("--cutoff", type=float, default=0, help="Minimum value cutoff for stacking")
    subparser.add_argument("--cols_ord", nargs="+", help="Color column value order")
    subparser.add_argument("--x_ord", nargs="+", help="X-axis value order")

    subparser.add_argument("--file", type=str, help="Output filename")
    subparser.add_argument("--dir", type=str, help="Output directory path")
    subparser.add_argument("--cmap", type=str, default="Set2", help="Colormap")

    subparser.add_argument("--errcap", type=int, default=4, help="Error bar cap width")
    subparser.add_argument("--vertical", action="store_true", help="Stack bars vertically (default True)")

    # Figure & layout
    subparser.add_argument("--figsize", nargs=2, type=int, default=(10, 6), help="Figure size (width height)")
    subparser.add_argument("--title", type=str, default="")
    subparser.add_argument("--title_size", type=int, default=18)
    subparser.add_argument("--title_weight", type=str, default="bold")
    subparser.add_argument("--title_font", type=str, default="Arial")

    # X-axis formatting
    subparser.add_argument("--x_axis", type=str, default="")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_axis_font", type=str, default="Arial")
    subparser.add_argument("--x_axis_dims", nargs=2, type=float, default=(0,0))
    subparser.add_argument("--x_ticks_rot", type=int, help="X-axis tick rotation")
    subparser.add_argument("--x_ticks_font", type=str, default="Arial")

    # Y-axis formatting
    subparser.add_argument("--y_axis", type=str, default="")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_axis_font", type=str, default="Arial")
    subparser.add_argument("--y_axis_dims", nargs=2, type=float, default=(0,0))
    subparser.add_argument("--y_ticks_rot", type=int, help="Y-axis tick rotation")
    subparser.add_argument("--y_ticks_font", type=str, default="Arial")

    # Legend options
    subparser.add_argument("--legend_title", type=str, default="")
    subparser.add_argument("--legend_title_size", type=int, default=12)
    subparser.add_argument("--legend_size", type=int, default=12)
    subparser.add_argument("--legend_bbox_to_anchor", nargs=2, type=float, default=(1, 1))
    subparser.add_argument("--legend_loc", type=str, default="upper left")
    subparser.add_argument("--legend_ncol", type=int, default=1)

    # Display and formatting
    subparser.add_argument("--show", action="store_true", help="Show the plot")
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize labels/legends with spacing")

'''
    add_common_plot_vol_args(subparser): Add common arguments for volcano plot
'''
def add_common_plot_vol_args(subparser):
    # Required arguments
    subparser.add_argument("--df", type=str, help="Input file")
    subparser.add_argument("--x", type=str, help="X-axis column (e.g. FC)")
    subparser.add_argument("--y", type=str, help="Y-axis column (e.g. pval)")

    # Optional data columns
    subparser.add_argument("--stys", type=str, help="Style column name")
    subparser.add_argument("--size", type=str, help="Size column name")
    subparser.add_argument("--size_dims", nargs=2, type=float, help="Min/max for size scaling")
    subparser.add_argument("--label", type=str, help="Label column name")

    # Thresholds
    subparser.add_argument("--FC_threshold", type=float, default=2, help="Fold change threshold")
    subparser.add_argument("--pval_threshold", type=float, default=0.05, help="P-value threshold")

    # Output
    subparser.add_argument("--file", type=str, help="Output file name")
    subparser.add_argument("--dir", type=str, help="Output directory path")

    # Aesthetics
    subparser.add_argument("--color", type=str, default="lightgray", help="Color for nonsignificant values")
    subparser.add_argument("--alpha", type=float, default=0.5, help="Transparency for nonsignificant values")
    subparser.add_argument("--edgecol", type=str, default="black", help="Edge color")
    subparser.add_argument("--vertical", action="store_true", help="Use vertical layout (default: True)")

    # Figure setup
    subparser.add_argument("--figsize", nargs=2, type=int, default=(10, 6), help="Figure size (width height)")
    subparser.add_argument("--title", type=str, default="")
    subparser.add_argument("--title_size", type=int, default=18)
    subparser.add_argument("--title_weight", type=str, default="bold")
    subparser.add_argument("--title_font", type=str, default="Arial")

    # X-axis settings
    subparser.add_argument("--x_axis", type=str, default="")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_axis_font", type=str, default="Arial")
    subparser.add_argument("--x_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--x_ticks_rot", type=int, default=0)
    subparser.add_argument("--x_ticks_font", type=str, default="Arial")
    subparser.add_argument("--x_ticks", nargs="+", help="Custom x-axis tick values")

    # Y-axis settings
    subparser.add_argument("--y_axis", type=str, default="")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_axis_font", type=str, default="Arial")
    subparser.add_argument("--y_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--y_ticks_rot", type=int, default=0)
    subparser.add_argument("--y_ticks_font", type=str, default="Arial")
    subparser.add_argument("--y_ticks", nargs="+", help="Custom y-axis tick values")

    # Legend
    subparser.add_argument("--legend_title", type=str, default="")
    subparser.add_argument("--legend_title_size", type=int, default=12)
    subparser.add_argument("--legend_size", type=int, default=9)
    subparser.add_argument("--legend_bbox_to_anchor", nargs=2, type=float, default=(1, 1))
    subparser.add_argument("--legend_loc", type=str, default="upper left")
    subparser.add_argument("--legend_items", nargs=2, type=int, default=(0, 0))
    subparser.add_argument("--legend_ncol", type=int, default=1)

    # Boolean switches
    subparser.add_argument("--display_size", action="store_true", help="Display size annotations")
    subparser.add_argument("--display_labels", action="store_true", help="Display text labels")
    subparser.add_argument("--return_df", action="store_true", help="Return annotated DataFrame")
    subparser.add_argument("--show", action="store_true", help="Show the plot")
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize labels/legends with spaces")

# Main method
'''
    main(): Endogenous Deep Mutational Scans
'''
def main():
    print("project: Endogenous Deep Mutational Scans")

    # Add parser and subparsers
    parser = argparse.ArgumentParser(description="Endogenous Deep Mutational Scans", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()

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
    parser_stat_difference.add_argument("--para", action="store_true", help="Use parametric test (default: True)")
    parser_stat_difference.add_argument("--alpha", type=float, default=0.05, help="Significance level (default: 0.05)")
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
                                         help="Correlation method to use (default: pearson)")
    parser_stat_correlation.add_argument("--numeric_only", action="store_true", help="Only use numeric columns (default: True)")

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
    '''
    parser_io = subparsers.add_parser("io", help="Input/Output")
    subparsers_io = parser_io.add_subparsers()
    
    # Create subparsers for commands
    parser_io_in_subs = subparsers_io.add_parser("in_subs", help="Moves all files with a given suffix into subfolders named after the files (excluding the suffix)")
    parser_io_out_subs = subparsers_io.add_parser("out_subs", help="Delete subdirectories and move their files to the parent directory")
    parser_io_create_sh = subparsers_io.add_parser("create_sh", help='Generate SLURM shell script for Harvard FASRC cluster.')
    parser_io_split_R1_R2 = subparsers_io.add_parser("split_R1_R2", help='Split paired reads into new R1 and R2 subdirectories at the parent directory.')
    
    # Add common arguments
    for parser_io_common in [parser_io_in_subs,parser_io_out_subs,parser_io_split_R1_R2]:
        parser_io_common.add_argument("--dir", help="Path to parent directory", type=str, required=True)
    
    # in_subs() arguments
    parser_io_in_subs.add_argument("--suf", help="File suffix (e.g., '.txt', '.csv') to filter files.", type=int, required=True) 
    
    # create_sh(): creates a shell script with SLURM job submission parameters for Harvard FASRC cluster.
    parser_io_create_sh.add_argument('--dir', type=str, required=True, help='Directory to save the shell script.')
    parser_io_create_sh.add_argument('--file', type=str, required=True, help='Name of the shell script file to create.')
    parser_io_create_sh.add_argument('--cores', type=int, default=1, help='Number of CPU cores to request.')
    parser_io_create_sh.add_argument('--partition', type=str, default='serial_requeue', help='SLURM partition to use.')
    parser_io_create_sh.add_argument('--time', type=str, default='0-00:10', help='Job run time in D-HH:MM format.')
    parser_io_create_sh.add_argument('--mem', type=int, default=1000, help='Memory in MB.')
    parser_io_create_sh.add_argument('--email', type=str, default=None, help='Notification email address.')
    parser_io_create_sh.add_argument('--python', type=str, default='python/3.12.5-fasrc01', help='Python module to load.')
    parser_io_create_sh.add_argument('--env', type=str, default='edms', help='Conda environment to activate.')

    # Call command functions
    parser_io_in_subs.set_defaults(func=io.in_subs)
    parser_io_out_subs.set_defaults(func=io.out_subs)
    parser_io_create_sh.set_defaults(func=io.create_sh)
    parser_io_split_R1_R2.set_defaults(func=io.split_R1_R2)

    '''
    edms.gen.cli:
    - access(): make all files and subdirectories accessible on Harvard FASRC
    - smaller_fastq(): create new subdirectory containing fastqs with the # of reads limited
    '''
    parser_cli = subparsers.add_parser("cli", help="Command Line Interaction")
    subparsers_cli = parser_cli.add_subparsers()
    
    # Create subparsers for commands
    parser_cli_access = subparsers_cli.add_parser("access", help="Make all files and subdirectories accessible on Harvard FASRC")
    parser_cli_smaller_fastq = subparsers_cli.add_parser("smaller_fastq", help="Create new subdirectory containing fastqs with the # of reads limited")
    
    # Add common arguments
    for parser_cli_common in [parser_cli_access, parser_cli_smaller_fastq]:
        parser_cli_common.add_argument("--pt", help="Path to parent directory", type=str, default='.')
    
    # Smaller_fastq arguments
    parser_cli_smaller_fastq.add_argument("--reads", help="# of reads per fastq file", type=int, default='100000') 
    parser_cli_smaller_fastq.add_argument("--suf", help="Fastq file suffix", type=int, default=".fastq.gz") 
    
    # Call command functions
    parser_cli_access.set_defaults(func=cli.access)
    parser_cli_smaller_fastq.set_defaults(func=cli.smaller_fastq)

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
    parser_ngs_pcrs.add_argument("--cycles", help="Number of cycles for PCR1", type=str, default='30')
    parser_ngs_pcrs.add_argument("--ultra", help="Using NEB Ultra II reagents", action="store_true")
    parser_ngs_pcrs.add_argument('--total_uL', type=int, default=20, help='Total reaction volume (uL)')
    parser_ngs_pcrs.add_argument('--mm_x', type=float, default=1.1, help='Master mix multiplier')
    
    # pcrs(): Column names
    parser_ngs_pcrs.add_argument('--gDNA_id_col', default='ID', help='gDNA ID column name')
    parser_ngs_pcrs.add_argument('--pcr1_id_col', default='PCR1 ID', help='PCR1 ID column name')
    parser_ngs_pcrs.add_argument('--pcr1_fwd_col', default='PCR1 FWD', help='PCR1 FWD column name')
    parser_ngs_pcrs.add_argument('--pcr1_rev_col', default='PCR1 REV', help='PCR1 REV column name')
    parser_ngs_pcrs.add_argument('--pcr2_id_col', default='PCR2 ID', help='PCR2 ID column name')
    parser_ngs_pcrs.add_argument('--pcr2_fwd_col', default='PCR2 FWD', help='PCR2 FWD column name')
    parser_ngs_pcrs.add_argument('--pcr2_rev_col', default='PCR2 REV', help='PCR2 REV column name')

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
    
    # hamming_distance_matrix(): compute pairwise Hamming distance matrix for a list of sequences stored in a dataframe
    parser_ngs_hamming = subparsers_ngs.add_parser("hamming", help="Compute pairwise Hamming distance matrix")
    
    parser_ngs_hamming.add_argument("--df", help="Input file", type=str, required=True)
    parser_ngs_hamming.add_argument("--id", help="ID column name", type=str, required=True)
    parser_ngs_hamming.add_argument("--seqs", help="Sequences column name", type=str, required=True)
    
    parser_ngs_hamming.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_ngs_hamming.add_argument("--file", help="Output file name", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_hamming.csv')
    
    parser_ngs_hamming.set_defaults(func=ngs.hamming_distance_matrix)

    '''
    edms.bio.clone
    - sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
    - epegRNAs(): design GG cloning oligonucleotides for prime editing epegRNAs
    - ngRNAs(): design GG cloning oligonucleotides for prime editing ngRNAs
    - ng_epegRNAs(): design GG cloning oligonucleotides for prime editing ng/epegRNAs (all-in-one vector)
    - pe_twist_oligos(): makes twist oligonucleotides for prime editing
    - pcr_sim(): returns dataframe with simulated pcr product 
    '''
    parser_clone = subparsers.add_parser("clone", help="Molecular cloning")
    subparsers_clone = parser_clone.add_subparsers()

    # sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
    parser_clone_sgRNAs = subparsers_clone.add_parser("sgRNAs", help="Design GG oligos for sgRNAs (cutting or BE)")
    
    parser_clone_sgRNAs.add_argument("--df", type=str, help="Input file path",required=True)
    parser_clone_sgRNAs.add_argument("--id", type=str, help="Column name for unique sgRNA identifier",required=True)

    parser_clone_sgRNAs.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_clone_sgRNAs.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_sgRNAs.csv')

    parser_clone_sgRNAs.add_argument("--tG", action="store_true", help="Add 5' G to spacer if needed")
    parser_clone_sgRNAs.add_argument("--order", action="store_true", help="Format output for ordering oligos")
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

    # ng_epegRNAs(): design GG cloning oligonucleotides for prime editing ng/epegRNAs (all-in-one vector)
    parser_clone_ng_epegRNAs = subparsers_clone.add_parser("ng_epegRNAs", help="Design GG oligos for ng/epegRNAs (all-in-one vector)")

    parser_clone_ng_epegRNAs.add_argument("--df", type=str, help="Input file path", required=True)
    parser_clone_ng_epegRNAs.add_argument("--id", type=str, help="Column name for unique sequence identifier",required=True)

    parser_clone_ng_epegRNAs.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_clone_ng_epegRNAs.add_argument("--file", help="Output file name", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_ng_epegRNAs.csv')

    parser_clone_ng_epegRNAs.add_argument("--dont_tG", dest="tG", default=True, action="store_false", help="Don't add 5' G to spacer if needed")
    parser_clone_ng_epegRNAs.add_argument("--dont_order", dest="order", default=True, action="store_false", help="Don't format output for ordering oligos")
    parser_clone_ng_epegRNAs.add_argument("--order_epegRNA_scaffold", action="store_true", help="Include epegRNA scaffold sequence in the oligo order")
    parser_clone_ng_epegRNAs.add_argument("--dont_make_extension", dest="make_extension", default=True, action="store_false", help="Don't build extension from RTT, PBS, and linker")
    parser_clone_ng_epegRNAs.add_argument("--ng_spacer", type=str, default="Spacer_sequence_ngRNA", help="Column name for ngRNA spacer sequence")
    parser_clone_ng_epegRNAs.add_argument("--ng_spacer_t5", type=str, default="GTTT", help="Top 5' overhang for ngRNA spacer")
    parser_clone_ng_epegRNAs.add_argument("--ng_spacer_t3", type=str, default="", help="Top 3' overhang for ngRNA spacer")
    parser_clone_ng_epegRNAs.add_argument("--ng_spacer_b5", type=str, default="AAAC", help="Bottom 5' overhang for ngRNA spacer")
    parser_clone_ng_epegRNAs.add_argument("--ng_spacer_b3", type=str, default="", help="Bottom 3' overhang for ngRNA spacer")
    parser_clone_ng_epegRNAs.add_argument("--epeg_spacer", type=str, default="Spacer_sequence_epegRNA", help="Column name for epegRNA spacer sequence")
    parser_clone_ng_epegRNAs.add_argument("--epeg_spacer_t5", type=str, default="CACC", help="Top 5' overhang for epegRNA spacer")
    parser_clone_ng_epegRNAs.add_argument("--epeg_spacer_t3", type=str, default="GTTTAAGAGC", help="Top 3' overhang for epegRNA spacer")
    parser_clone_ng_epegRNAs.add_argument("--epeg_spacer_b5", type=str, default="", help="Bottom 5' overhang for epegRNA spacer")
    parser_clone_ng_epegRNAs.add_argument("--epeg_spacer_b3", type=str, default="", help="Bottom 3' overhang for epegRNA spacer")
    parser_clone_ng_epegRNAs.add_argument("--extension", type=str, default="Extension_sequence", help="Column name for epegRNA extension sequence")
    parser_clone_ng_epegRNAs.add_argument("--extension_t5", type=str, default="", help="Top 5' overhang for epegRNA extension")
    parser_clone_ng_epegRNAs.add_argument("--extension_t3", type=str, default="", help="Top 3' overhang for epegRNA extension")
    parser_clone_ng_epegRNAs.add_argument("--extension_b5", type=str, default="CGCG", help="Bottom 5' overhang for epegRNA extension")
    parser_clone_ng_epegRNAs.add_argument("--extension_b3", type=str, default="GCACCGACTC", help="Bottom 3' overhang for epegRNA extension")
    parser_clone_ng_epegRNAs.add_argument("--RTT", type=str, default="RTT_sequence", help="Column name for RTT (reverse transcriptase template)")
    parser_clone_ng_epegRNAs.add_argument("--PBS", type=str, default="PBS_sequence", help="Column name for PBS (primer binding site)")
    parser_clone_ng_epegRNAs.add_argument("--linker", type=str, default="Linker_sequence", help="Column name for linker")

    parser_clone_ng_epegRNAs.set_defaults(func=cl.ng_epegRNAs)

    # pe_twist_oligos(): makes twist oligonucleotides for prime editing
    parser_clone_pe_twist = subparsers_clone.add_parser("pe_twist", help="Design Twist oligos for PE constructs")
    
    parser_clone_pe_twist.add_argument("--df", type=str, help="Input file path", required=True)
    parser_clone_pe_twist.add_argument("--id_pre", type=str, help="Prefix for ID column", required=True)

    parser_clone_pe_twist.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_clone_pe_twist.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pe_twist.csv')

    parser_clone_pe_twist.add_argument("--tG", action="store_true", help="Add 5' G to spacers if needed")
    parser_clone_pe_twist.add_argument("--make_extension", action="store_true", help="Build extension from RTT, PBS, and linker")
    parser_clone_pe_twist.add_argument("--UMI_length", type=int, default=8, help="Length of UMI (default: 8)")
    parser_clone_pe_twist.add_argument("--UMI_GC_fract", nargs=2, type=float, default=(0.4, 0.6), help="Tuple for GC content bounds (e.g. 0.4 0.6)")
    parser_clone_pe_twist.add_argument("--fwd_barcode_t5", type=str, default="Forward Barcode", help="Forward barcode column name")
    parser_clone_pe_twist.add_argument("--rev_barcode_t3", type=str, default="Reverse Barcode", help="Reverse barcode column name")
    parser_clone_pe_twist.add_argument("--homology_arm_t5", type=str, default="Homology Arm 5", help="Homology arm 5' column name")
    parser_clone_pe_twist.add_argument("--homology_arm_t3", type=str, default="Homology Arm 3", help="Homology arm 3' column name")
    parser_clone_pe_twist.add_argument("--ngRNA_hU6_gg_insert", type=str, default="GTTTAGAGACGATCGACGTCTCACACC", help="Insert sequence for hU6 Golden Gate ngRNA")
    parser_clone_pe_twist.add_argument("--epegRNA_gg_insert", type=str, default="GTTTAAGAGCAGGTGCTAGACCTGCGTCGGTGC", help="Insert sequence for Golden Gate epegRNA")
    parser_clone_pe_twist.add_argument("--ngRNA_spacer", type=str, default="Spacer_sequence_ngRNA", help="ngRNA spacer column")
    parser_clone_pe_twist.add_argument("--epegRNA_spacer", type=str, default="Spacer_sequence_epegRNA", help="epegRNA spacer column")
    parser_clone_pe_twist.add_argument("--epegRNA_extension", type=str, default="Extension_sequence", help="epegRNA extension column")
    parser_clone_pe_twist.add_argument("--epegRNA_RTT", type=str, default="RTT_sequence", help="RTT column name")
    parser_clone_pe_twist.add_argument("--epegRNA_PBS", type=str, default="PBS_sequence", help="PBS column name")
    parser_clone_pe_twist.add_argument("--epegRNA_linker", type=str, default="Linker_sequence", help="Linker column name")
    parser_clone_pe_twist.add_argument("--epegRNA_pbs_length", type=str, default="PBS_length", help="PBS length column name")
    parser_clone_pe_twist.add_argument("--ngRNA_group", type=str, default="ngRNA_group", help="Group column name for ngRNAs")

    parser_clone_pe_twist.set_defaults(func=cl.pe_twist_oligos) 
    
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
    '''
    parser_fastq = subparsers.add_parser("fastq", help="FASTQ files")
    subparsers_fastq = parser_fastq.add_subparsers()

    parser_fastq_revcom = subparsers_fastq.add_parser("revcom", help="Reverse complement all FASTQ files in a directory")
    parser_fastq_unzip = subparsers_fastq.add_parser("unzip", help="Unzip gzipped FASTQ files to a new directory")
    parser_fastq_comb = subparsers_fastq.add_parser("comb", help="Combine multiple FASTQ files into a single FASTQ (.fastq or .fastq.gz)")
    parser_fastq_genotyping = subparsers_fastq.add_parser("genotyping", help="Quantify edit outcomes from sequencing data")

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
    
    parser_fastq_genotyping.add_argument("--desired_edits", type=str, nargs='+', help="List of desired edits", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--qall", type=int, help="Minimum Phred quality score for all bases", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--qtrim", type=int, help="Phred quality threshold for end trimming", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--qavg", type=int, help="Minimum average Phred quality score", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--qmask", type=int, help="Phred quality threshold for masking to N", default=argparse.SUPPRESS)

    parser_fastq_genotyping.add_argument("--no_saves", action="store_false", help="Don't save read statistics and genotypes files",dest="save",default=True)
    parser_fastq_genotyping.add_argument("--masks", action="store_true", help="Include masked sequence and translation",default=False)
    parser_fastq_genotyping.add_argument("--keepX", action="store_true", help="Keep unknown translation (X) in output", default=False)
    
    parser_fastq_revcom.set_defaults(func=fq.revcom_fastqs)
    parser_fastq_unzip.set_defaults(func=fq.unzip_fastqs)
    parser_fastq_comb.set_defaults(func=fq.comb_fastqs)
    parser_fastq_genotyping.set_defaults(func=fq.genotyping)

    '''
    Add pe.py
    - epegRNA_linkers(): generate epegRNA linkers between PBS and 3' hairpin motif & finish annotations
    
    WIP: Full workflow (DO THIS IF PRIME DESIGN WORKS BETTER THAN OptiPrime)
    1. PrimeDesigner(): 
        - PrimeDesignInput(): creates and checks PrimeDesign saturation mutagenesis input file
        - PrimeDesign(): run PrimeDesign using Docker (NEED TO BE RUNNING DESKTOP APP)
        - PrimeDesignOutput(): splits peg/ngRNAs from PrimeDesign output & finishes annotations
    2. Pilot_Screen():
        - shared_sequences(): Reduce PE library into shared spacers and PBS sequences
        - cosmic/cvar: priority_muts(), priority_edits()
    3. epegRNA_linkers(): generate epegRNA linkers between PBS and 3' hairpin motif & finish annotations
    4. MergePrimeDesignOutput(): rejoins epeg/ngRNAs from PrimeDesign output & creates ngRNA_groups
    '''
    parser_pe = subparsers.add_parser("pe", help="Prime Editing")
    subparsers_pe = parser_pe.add_subparsers()

    parser_pe_epegRNA_linkers = subparsers_pe.add_parser("epegRNA_linkers", help="Generate epegRNA linkers between PBS and 3' hairpin motif")
    
    # epegRNA_linkers(): Required input file
    parser_pe_epegRNA_linkers.add_argument('--pegRNAs', help='Path to pegRNAs file',required=True)

    # epegRNA_linkers(): Optional parameters
    parser_pe_epegRNA_linkers.add_argument('--epegRNA_motif_sequence', default='CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA', help='epegRNA motif sequence (default: tevopreQ1)')
    parser_pe_epegRNA_linkers.add_argument('--checkpoint_dir', type=str, help='Checkpoint directory path')
    parser_pe_epegRNA_linkers.add_argument('--checkpoint_file', help='Checkpoint file name')
    parser_pe_epegRNA_linkers.add_argument('--checkpoint_pt', type=str, default='', help='Previous checkpoint full path')
    
    # Set defaults
    parser_pe_epegRNA_linkers.set_defaults(func=pe.epegRNA_linkers)

    # Parse all arguments
    args = parser.parse_args()
    args_dict = vars(args)
    func = args_dict.pop("func")
    func(**args_dict)    