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
from .gen import plot as p
from .gen import com
from .bio import ngs

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

    # X-axis settings
    subparser.add_argument("--x_axis", type=str, default="", help="X-axis label")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_axis_scale", type=str, default="linear")
    subparser.add_argument("--x_axis_dims", nargs=2, type=float, default=(0,0))
    subparser.add_argument("--x_ticks_rot", type=int, default=0)
    subparser.add_argument("--xticks", nargs="+", help="X-axis ticks")

    # Y-axis settings
    subparser.add_argument("--y_axis", type=str, default="", help="Y-axis label")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_axis_scale", type=str, default="linear")
    subparser.add_argument("--y_axis_dims", nargs=2, type=float, default=(0,0))
    subparser.add_argument("--y_ticks_rot", type=int, default=0)
    subparser.add_argument("--yticks", nargs="+", help="Y-axis ticks")

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

    # X-axis settings
    subparser.add_argument("--x_axis", type=str, default="")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_axis_scale", type=str, default="linear")
    subparser.add_argument("--x_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--x_ticks_rot", type=int, default=0)
    subparser.add_argument("--xticks", nargs="+", help="X-axis ticks")

    # Y-axis settings
    subparser.add_argument("--y_axis", type=str, default="")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_axis_scale", type=str, default="linear")
    subparser.add_argument("--y_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--y_ticks_rot", type=int, default=0)
    subparser.add_argument("--yticks", nargs="+", help="Y-axis ticks")

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

    # X-axis
    subparser.add_argument("--x_axis", type=str, default="")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_axis_scale", type=str, default="linear")
    subparser.add_argument("--x_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--x_ticks_rot", type=int, default=0)
    subparser.add_argument("--xticks", nargs="+", help="X-axis tick values")

    # Y-axis
    subparser.add_argument("--y_axis", type=str, default="")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_axis_scale", type=str, default="linear")
    subparser.add_argument("--y_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--y_ticks_rot", type=int, default=0)
    subparser.add_argument("--yticks", nargs="+", help="Y-axis tick values")

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
    subparser.add_argument("--x", type=str, help="X-axis column", required=True)
    subparser.add_argument("--y", type=str, help="Y-axis column", required=True)

    # Optional arguments
    subparser.add_argument("--vars", type=str, default="variable", help="Name of variable column")
    subparser.add_argument("--vals", type=str, default="value", help="Name of value column")
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
    subparser.add_argument("--figsize", nargs=2, type=int, default=(10, 6), help="Figure size (width height)")

    # X-axis
    subparser.add_argument("--x_axis", type=str, default="")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_ticks_rot", type=int, default=0)

    # Y-axis
    subparser.add_argument("--y_axis", type=str, default="")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_ticks_rot", type=int, default=0)

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

    # X-axis formatting
    subparser.add_argument("--x_axis", type=str, default="")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_ticks_rot", type=int, help="X-axis tick rotation")
    subparser.add_argument("--x_ticks_ha", type=str, help="X-axis tick horizontal alignment")

    # Y-axis formatting
    subparser.add_argument("--y_axis", type=str, default="")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_ticks_rot", type=int, help="Y-axis tick rotation")
    subparser.add_argument("--y_ticks_ha", type=str, help="Y-axis tick horizontal alignment")
    subparser.add_argument("--yticks", nargs="+", help="Custom y-axis tick values")

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

    # X-axis settings
    subparser.add_argument("--x_axis", type=str, default="")
    subparser.add_argument("--x_axis_size", type=int, default=12)
    subparser.add_argument("--x_axis_weight", type=str, default="bold")
    subparser.add_argument("--x_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--x_ticks_rot", type=int, default=0)
    subparser.add_argument("--xticks", nargs="+", help="Custom x-axis tick values")

    # Y-axis settings
    subparser.add_argument("--y_axis", type=str, default="")
    subparser.add_argument("--y_axis_size", type=int, default=12)
    subparser.add_argument("--y_axis_weight", type=str, default="bold")
    subparser.add_argument("--y_axis_dims", nargs=2, type=float, default=(0, 0))
    subparser.add_argument("--y_ticks_rot", type=int, default=0)
    subparser.add_argument("--yticks", nargs="+", help="Custom y-axis tick values")

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
    NEED TO ADD IO.SAVE CONDITIONAL WITH DIR AND FILE
    '''



    '''
    edms.gen.com:
    - expand_subs(): delete subdirectories and move their files to the parent directory
    - split_paired_reads(): split paired reads into new R1 and R2 subdirectories at the parent directory
    - smaller_fastq(): create new subdirectory containing fastqs with the # of reads limited
    '''
    parser_com = subparsers.add_parser("com", help="Command Line Interaction")
    subparsers_com = parser_com.add_subparsers()
    
    # Create subparsers for commands
    parser_com_expand_subs = subparsers_com.add_parser("expand_subs", help="Delete subdirectories and move their files to the parent directory")
    parser_com_split_paired_reads = subparsers_com.add_parser("split_paired_reads", help="Split paired reads into new R1 and R2 subdirectories at the parent directory")
    parser_com_smaller_fastq = subparsers_com.add_parser("smaller_fastq", help="Ccreate new subdirectory containing fastqs with the # of reads limited")
    
    # Add common arguments
    for parser_com_common in [parser_com_expand_subs, parser_com_split_paired_reads, parser_com_smaller_fastq]:
        parser_com_common.add_argument("--pt", help="Path to parent directory", type=str, required=True)
    
    # Smaller_fastq arguments
    parser_com_smaller_fastq.add_argument("--reads", help="# of reads per fastq file", type=int, required=True) 
    parser_com_smaller_fastq.add_argument("--suf", help="Fastq file suffix", type=int, default=".fastq.gz") 
    
    # Call command functions
    parser_com_expand_subs.set_defaults(func=com.expand_subs)
    parser_com_split_paired_reads.set_defaults(func=com.split_paired_reads)
    parser_com_smaller_fastq.set_defaults(func=com.smaller_fastq)

    '''
    edms.dat.cosmic:
    '''

    '''
    edms.dat.cvar:
    '''

    '''
    edms.bio.ngs:
    - pcrs(): generates NGS PCR plan automatically (Default: 96-well plates including outer wells)
    - compute_distance_matrix(): compute pairwise Hamming distance matrix for a list of sequences stored in a dataframe
    '''

    parser_ngs = subparsers.add_parser("ngs", help="Next generation sequencing")
    subparsers_ngs = parser_ngs.add_subparsers()

    # pcrs(): generates NGS PCR plan automatically (Default: 96-well plates including outer wells)
    parser_ngs_pcrs = subparsers_ngs.add_parser("pcrs", help="Plan NGS PCRs")
    parser_ngs_pcrs.add_argument("--df", help="Input file", type=str, required=True)
    parser_ngs_pcrs.add_argument("--dir", help="Output directory path", type=str, default='.')
    parser_ngs_pcrs.add_argument("--file", help="Output file name (.xlsx)", type=str, default='NGS_plan.xlsx')
    parser_ngs_pcrs.set_defaults(func=ngs.pcrs)
    
    # hamming_distance_matrix(): compute pairwise Hamming distance matrix for a list of sequences stored in a dataframe
    parser_ngs_hamming = subparsers_ngs.add_parser("hamming", help="Compute pairwise Hamming distance matrix")
    parser_ngs_hamming.add_argument("--df", help="Input file", type=str, required=True)
    parser_ngs_hamming.add_argument("--id", help="ID column name", type=str, required=True)
    parser_ngs_hamming.add_argument("--seqs", help="Sequences column name", type=str, required=True)
    parser_ngs_hamming.add_argument("--dir", help="Output directory path", type=str, default='.')
    parser_ngs_hamming.add_argument("--file", help="Output file name", type=str, default='hamming.csv')
    
    '''
    edms.bio.clone
    '''

    '''
    edms.bio.transfect
    '''

    '''
    edms.bio.qPCR:
    '''
    
    '''
    edms.bio.fastq:
    '''

    # Parse all arguments
    args = parser.parse_args()
    args_dict = vars(args)
    func = args_dict.pop("func")
    func(**args_dict)    