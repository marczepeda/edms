''' 
src/edms/gen/cli.py             Command Line Interface for EDMS General module
├── __init__.py                 Initializer
├── html.py                     HTML module
├── io.py                       Input/Output module
├── plot.py                     Plot module
├── stat.py                     Statistical Analysis module
├── tidy.py                     Tidy Data module
├── com.py                      Command-line Interaction module
└── image.py                    Image Processing module

Usage:
[fastq helper methods]
- add_common_fastq_label_args(subparser): Add common label arguments for fastq-related subparsers

[Plot subparser methods]
- add_common_plot_scat_args(subparser, fastq_parser=False): Add common arguments for scatter plot related graphs
- add_common_plot_cat_args(subparser): Add common arguments for categorical graphs
- add_common_plot_dist_args(subparser): Add common arguments for distribution graphs
- add_common_plot_heat_args(subparser): Add common arguments for heatmap graphs
- add_common_plot_stack_args(subparser): Add common arguments for stacked bar plot
- add_common_plot_vol_args(subparser, fastq_parser=False): Add common arguments for volcano plot

[Main subparser method]
- add_subparser(): Attach all gen-related subparsers to the top-level CLI.
'''
import argparse
import datetime
import sys # might use later
from rich import print as rprint # might use later

from . import com, io, plot as p, stat as st, html as ht
from ..utils import parse_tuple_int, parse_tuple_float # might use later

# fastq helper methods
def add_common_fastq_label_args(subparser):
    '''
    add_common_fastq_label_args(subparser): Add common label arguments for fastq-related subparsers
    '''
    subparser.add_argument("--label_size", type=int, help="Label font size (Default: 16)", default=16)
    subparser.add_argument("--no_label_info", dest='label_info', action="store_false", help="Don't include additional info for labels if .html plot (Default: True)", default=True)
    subparser.add_argument("--aa_properties", nargs="+", help="Use aa_properties to format labels (Options: hydrophobicity, polarity, charge, vdw_volume, pKa_C_term, pKa_N_term, pKa_side_chain)", default=argparse.SUPPRESS)
    subparser.add_argument("--cBioPortal", type=str, help="Gene name (if saved to ~/.config/edms/cBioPortal_mutations) or file path for cBioPortal mutation data processed through edms.dat.cBioPortal.mutations()", default=argparse.SUPPRESS)
    subparser.add_argument("--UniProt", type=str, help="UniProt accession (if saved to ~/.config/edms/UniProt) or file path for UniProt flat file. See edms.dat.uniprot.retrieve_flat_file() or edms uniprot retrieve -h for more information.", default=argparse.SUPPRESS)
    subparser.add_argument("--PhosphoSitePlus", type=str, help="UniProt accession", default=argparse.SUPPRESS)
    subparser.add_argument("--PDB_contacts", type=str, help="PDB ID (if saved to ~/.config/edms/PDB) or file path for PDB structure file. See edms.dat.pdb.retrieve() or edms uniprot retrieve -h for more information.", default=argparse.SUPPRESS)
    subparser.add_argument("--PDB_neighbors", type=str, help="PDB ID (if saved to ~/.config/edms/PDB) or file path for PDB structure file. See edms.dat.pdb.retrieve() or edms uniprot retrieve -h for more information.", default=argparse.SUPPRESS)

# Plot subparser methods
def add_common_plot_scat_args(subparser, fastq_torn_parser=False, fastq_corr_parser=False):
    '''
    add_common_plot_scat_args(subparser): Add common arguments for scatter plot related graphs
    '''
    # scat(): Required arguments
    subparser.add_argument("--df", help="Input dataframe file path", type=str, required=True)
    if fastq_torn_parser == False and fastq_corr_parser == False:
        subparser.add_argument("--x", help="X-axis column", type=str, required=True)
        subparser.add_argument("--y", help="Y-axis column", type=str, required=True)
    else:
        if fastq_corr_parser == True:
            subparser.add_argument("--cond_col", help="condition column name for comparison", type=str, required=True)
            subparser.add_argument("--cond_vals", nargs="+", help="two condition values for comparison (x and y-axis)", type=str, required=True)
        subparser.add_argument("--FC", help="Fold change column name (Y-axis)", type=str, required=True)
        subparser.add_argument("--pval", help="p-value column name (size column if not specified)", type=str, required=True)
        
    # Optional core arguments
    if fastq_torn_parser == False and fastq_corr_parser == False:
        subparser.add_argument("--cols", type=str, help="Color column name")
        subparser.add_argument("--cols_ord", nargs="+", help="Column order (list of values)")
        subparser.add_argument("--cols_exclude", nargs="+", help="Columns to exclude from coloring")
        subparser.add_argument("--stys", type=str, help="Style column name")
        subparser.add_argument("--stys_order", nargs="+", help="Style order (list of values)")
        subparser.add_argument("--mark_order", nargs="+", help="Marker order (list of marker styles)")
    else:
        subparser.add_argument("--size", type=str, help="Column name used to scale point sizes (Default: -log10('pval'); specify 'false' for no size)")
        subparser.add_argument("--size_dims", type=parse_tuple_float, help="Size range for points formatted as min,max")
        if fastq_corr_parser == True:
            subparser.add_argument("--method", help="Correlation method (Default: 'pearson')", choices=['pearson', 'spearman', 'kendall'], type=str, default='pearson')
            subparser.add_argument("--not_weighted", dest='weighted', action='store_false', help="Weighted correlation by size column (Default: True)", default=True)

    subparser.add_argument("--label", type=str, help="Column name for point labels; static text for images, interactive tooltips for HTML")

    # Additional annotation data sources
    if fastq_torn_parser==True or fastq_corr_parser==True:
        add_common_fastq_label_args(subparser)
        if fastq_torn_parser==True:
            subparser.add_argument("--ss_h", type=int, help="Height for secondary structure in the plot (Default: autogenerate)", default=argparse.SUPPRESS)
            subparser.add_argument("--ss_y", type=int, help="Y position for secondary structure track in the plot (Default: autogenerate)", default=argparse.SUPPRESS)

    subparser.add_argument("--dir", help="Output directory path", type=str, default='./out')
    subparser.add_argument("--file", help="Output file name", type=str, required=False, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_scat.png')
    if fastq_torn_parser == False and fastq_corr_parser == False:
        subparser.add_argument("--palette_or_cmap", type=str, default="colorblind", help="Seaborn palette or matplotlib colormap")
    subparser.add_argument("--edgecol", type=str, default="black", help="Edge color for scatter points")

    # Figure appearance
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(5,5), help="Figure size as a tuple: width,height")
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
    subparser.add_argument("--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("--x_ticks_size", type=int, default=9, help="X-axis tick labels font size")
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
    subparser.add_argument("--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("--y_ticks_size", type=int, default=9, help="Y-axis tick labels font size")
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
    subparser.add_argument('--legend_columnspacing', type=int, default=argparse.SUPPRESS, help='space between columns in legend; only for html plots')
    subparser.add_argument('--legend_handletextpad', type=float, default=argparse.SUPPRESS, help='space between marker and text in legend; only for html plots')
    subparser.add_argument('--legend_labelspacing', type=float, default=argparse.SUPPRESS, help='vertical space between entries in legend; only for html plots')
    subparser.add_argument('--legend_borderpad', type=float, default=argparse.SUPPRESS, help='padding inside legend box; only for html plots')
    subparser.add_argument('--legend_handlelength', type=float, default=argparse.SUPPRESS, help='marker length in legend; only for html plots')
    subparser.add_argument('--legend_size_html_multiplier', type=float, default=argparse.SUPPRESS, help='legend size multiplier for html plots')

    # Display and formatting
    subparser.add_argument("--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
    subparser.add_argument("--show", action="store_true", help="Show the plot", default=False)
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize label/legend strings and replace underscores with spaces")
    if fastq_torn_parser == True or fastq_corr_parser == True:
        subparser.add_argument("--dont_display_legend", dest='display_legend', action="store_false", default=True, help="Display legend on plot (Default: True)")
        subparser.add_argument("--dont_display_labels", dest='display_labels', action="store_false", default=True, help="Display labels for significant values (Default: True)")
        subparser.add_argument("--dont_display_axis", dest='display_axis', action="store_false", default=True, help="Display x- and y-axis lines (Default: True)")

def add_common_plot_cat_args(subparser, fastq_parser=False):
    '''
    add_common_plot_cat_args(subparser): Add common arguments for categorical graphs
    '''
    # cat(): Required arguments
    if fastq_parser == True:
        subparser.add_argument("--type", dest="typ", help="Type of category plot", type=str, required=True, choices=['bar', 'box', 'violin', 'strip', 'swarm', 'point', 'count', 'bar_strip', 'box_strip', 'violin_strip', 'bar_swarm', 'box_swarm', 'violin_swarm'])
    subparser.add_argument("--df", help="Input dataframe file path", type=str, required=True)

    # Optional core arguments
    subparser.add_argument("--x", help="X-axis column name", type=str, default="")
    subparser.add_argument("--y", help="Y-axis column name", type=str, default="")
    subparser.add_argument("--cats_ord", nargs="+", help="Category column values order (x- or y-axis)")
    subparser.add_argument("--cats_exclude", nargs="+", help="Category column values exclude (x- or y-axis)")
    subparser.add_argument("--cols", type=str, help="Color column name for grouping")
    subparser.add_argument("--cols_ord", nargs="+", help="Color column values order")
    subparser.add_argument("--cols_exclude", nargs="+", help="Color column values to exclude")
    if fastq_parser == True:
        subparser.add_argument("--PDB_pt", type=str, help="PDB ID (if saved to ~/.config/edms/PDB) or file path for PDB structure file. See edms.dat.pdb.retrieve() or edms uniprot retrieve -h for more information", default=argparse.SUPPRESS)

    subparser.add_argument("--dir", type=str, help="Output directory", default='./out')
    subparser.add_argument("--file", type=str, help="Output filename", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_cat.png')
    subparser.add_argument("--palette_or_cmap", type=str, default="colorblind", help="Seaborn color palette or matplotlib colormap")
    subparser.add_argument("--edgecol", type=str, default="black", help="Edge color for markers")

    # Error bar and style options
    subparser.add_argument("--lw", type=int, default=1, help="Line width for plot edges")
    subparser.add_argument("--errorbar", type=str, default="sd", help="Error bar type: sd (standard deviation), etc.")
    subparser.add_argument("--errwid", type=float, default=1, help="Width of the error bars")
    subparser.add_argument("--errcap", type=float, default=0.1, help="Cap size on error bars")

    # Figure appearance
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(5,5), help="Figure size formatted 'width,height'")
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
    subparser.add_argument("--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("--x_ticks_size", type=int, default=9, help="Font size for X-axis tick labels")
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
    subparser.add_argument("--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("--y_ticks_size", type=int, default=9, help="Font size for Y-axis tick labels")
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
    subparser.add_argument('--legend_columnspacing', type=int, default=argparse.SUPPRESS, help='space between columns in legend; only for html plots')
    subparser.add_argument('--legend_handletextpad', type=float, default=argparse.SUPPRESS, help='space between marker and text in legend; only for html plots')
    subparser.add_argument('--legend_labelspacing', type=float, default=argparse.SUPPRESS, help='vertical space between entries in legend; only for html plots')
    subparser.add_argument('--legend_borderpad', type=float, default=argparse.SUPPRESS, help='padding inside legend box; only for html plots')
    subparser.add_argument('--legend_handlelength', type=float, default=argparse.SUPPRESS, help='marker length in legend; only for html plots')
    subparser.add_argument('--legend_size_html_multiplier', type=float, default=argparse.SUPPRESS, help='legend size multiplier for html plots')

    # Display and formatting
    subparser.add_argument("--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
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
    subparser.add_argument("--despine", action="store_true", help="Remove plot spines (despine)", default=False)

    # Figure appearance
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(5,5), help="Figure size formatted as 'width,height'")
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
    subparser.add_argument("--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("--x_ticks_size", type=int, default=9, help="Font size for X-axis tick labels")
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
    subparser.add_argument("--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("--y_ticks_size", type=int, default=9, help="Font size for Y-axis tick labels")
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
    subparser.add_argument("--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
    subparser.add_argument("--show", action="store_true", help="Show the plot in an interactive window", default=False)
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize and space legend/label values", default=False)

def add_common_plot_heat_args(subparser, fastq_parser=False, stat_parser=False):
    '''
    add_common_plot_heat_args(subparser): Add common arguments for heatmap graphs
    '''
    # Required arguments
    if stat_parser == False:
        subparser.add_argument("--df", help="Input dataframe file path", type=str, required=True)

    # Optional arguments
    if fastq_parser == False and stat_parser == False:
        subparser.add_argument("--x", type=str, help="X-axis column name to pivot tidy-formatted dataframe into matrix format")
        subparser.add_argument("--y", type=str, help="Y-axis column name to pivot tidy-formatted dataframe into matrix format")
        subparser.add_argument("--vars", type=str, help="Variable column name to split tidy-formatted dataframe into a dictionary of pivoted dataframes")
        subparser.add_argument("--vals", type=str, help="Value column name to populate pivoted dataframes")
    if fastq_parser == False:    
        subparser.add_argument("--vals_dims", type=parse_tuple_float, help="Value column limits formatted as 'vmin,vmax'")
    if fastq_parser == True:
        subparser.add_argument("--cond_col", type=str, help="Condition column name", required=True)
        subparser.add_argument("--cond", type=str, help="Condition value for filtering", required=True)
        subparser.add_argument("--FC", type=str, help="Fold change column name (values within heatmap after log2 transformation)", required=True)
        subparser.add_argument("--wt_prot", type=str, help="WT protein sequence", required=True)
        subparser.add_argument("--wt_res", type=int, help="WT protein sequence residue start number", required=True)
        subparser.add_argument("--cutoff", type=float, help="Comparison count mean cutoff for masking low-abundance values", default=argparse.SUPPRESS)
        subparser.add_argument("--aa", type=str, help="AA saturation mutagenesis. The 'aa' option [default] makes all amino acid substitutions ('aa_subs'), +1 amino acid insertions ('aa_ins'), and -1 amino acid deletions ('aa_dels').", choices=['aa', 'aa_subs', 'aa_ins', 'aa_dels'], default='aa')
        subparser.add_argument("--label", type=str, help="Label column name (Default: 'Edit'). Can't be None.", default='Edit')
    
    if stat_parser == False:
        subparser.add_argument("--dir", type=str, help="Output directory path", default='./out')
        subparser.add_argument("--file", type=str, help="Output filename", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_heat.png')
    subparser.add_argument("--edgecol", type=str, default="black", help="Color of cell edges")
    subparser.add_argument("--lw", type=int, default=1, help="Line width for cell borders")

    if fastq_parser == False:
        subparser.add_argument("--annot", action="store_true", help="Display cell values as annotations", default=False)
    subparser.add_argument("--center", type=float, default=0, help="Center value for colormap (Default: 0)")
    subparser.add_argument("--cmap", type=str, default="Reds", help="Matplotlib colormap to use for heatmap")
    if fastq_parser == True:
        subparser.add_argument("--cmap_WT", type=str, default="forestgreen", help="Color for WT values in colorbar (Default: forestgreen)")
        subparser.add_argument("--cmap_not_WT", type=str, default="lightgray", help="Color for non-WT values in colorbar (Default: lightgray)")
    subparser.add_argument("--sq", action="store_true", help="Use square aspect ratio for cells", default=False)
    subparser.add_argument("--cbar", action="store_true", help="Display colorbar", default=False)
    if stat_parser == False:
        subparser.add_argument("--cbar_label", type=str, default=argparse.SUPPRESS, help="Colorbar label")
    subparser.add_argument("--cbar_label_size", type=int, default=argparse.SUPPRESS, help="Font size for colorbar label")
    subparser.add_argument("--cbar_label_weight", type=str, default="bold", help="Font weight for colorbar label (Default: bold)", choices=['bold', 'normal', 'heavy'])
    subparser.add_argument("--cbar_tick_size", type=int, default=argparse.SUPPRESS, help="Font size for colorbar ticks")
    subparser.add_argument("--cbar_shrink", type=float, default=argparse.SUPPRESS, help="Shrink factor for colorbar")
    subparser.add_argument("--cbar_aspect", type=int, default=argparse.SUPPRESS, help="Aspect ratio for colorbar")
    subparser.add_argument("--cbar_pad", type=float, default=argparse.SUPPRESS, help="Padding for colorbar")
    subparser.add_argument("--cbar_orientation", type=str, default=argparse.SUPPRESS, help="Orientation of colorbar (Default: 'vertical')", choices=['vertical', 'horizontal'])

    # Title and size
    subparser.add_argument("--title", type=str, default="", help="Plot title")
    subparser.add_argument("--title_size", type=int, default=18, help="Font size of the title")
    subparser.add_argument("--title_weight", type=str, default="bold", help="Font weight of the title (e.g., bold, normal)")
    subparser.add_argument("--title_font", type=str, default="Arial", help="Font family for the title")
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(5,5), help="Figure size formatted as 'width,height'")

    # X-axis
    subparser.add_argument("--x_axis", type=str, default="", help="X-axis label")
    subparser.add_argument("--x_axis_size", type=int, default=12, help="Font size for X-axis label")
    subparser.add_argument("--x_axis_weight", type=str, default="bold", help="Font weight for X-axis label")
    subparser.add_argument("--x_axis_font", type=str, default="Arial", help="Font family for X-axis label")
    subparser.add_argument("--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("--x_ticks_size", type=int, default=9, help="Font size for X-axis tick labels")
    subparser.add_argument("--x_ticks_rot", type=int, default=0, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")

    # Y-axis
    subparser.add_argument("--y_axis", type=str, default="", help="Y-axis label")
    subparser.add_argument("--y_axis_size", type=int, default=12, help="Font size for Y-axis label")
    subparser.add_argument("--y_axis_weight", type=str, default="bold", help="Font weight for Y-axis label")
    subparser.add_argument("--y_axis_font", type=str, default="Arial", help="Font family for Y-axis label")
    subparser.add_argument("--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("--y_ticks_size", type=int, default=9, help="Font size for Y-axis tick labels")
    subparser.add_argument("--y_ticks_rot", type=int, default=0, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")

    # Final display
    subparser.add_argument("--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
    subparser.add_argument("--show", action="store_true", help="Show the plot in an interactive window", default=False)
    if fastq_parser == False:
        subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize and space labels/legend values", default=False)

def add_common_plot_stack_args(subparser, fastq_parser=False):
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
    
    subparser.add_argument("--cutoff_group", type=str, default=argparse.SUPPRESS, help="Column name to group by when applying cutoff")
    subparser.add_argument("--cutoff_value", type=float, default=0, help="Y-axis values needs be greater than (e.g. 0)")
    subparser.add_argument("--cutoff_remove", dest="cutoff_keep",action="store_false", help="Remove values below cutoff", default=True)
    subparser.add_argument("--cols_ord", nargs="+", help="Order of values in the color column")
    subparser.add_argument("--x_ord", nargs="+", help="Custom order of X-axis categories")
    subparser.add_argument("--palette_or_cmap", type=str, default="Set2", help="Seaborn palette or Matplotlib colormap for stacked bars")
    subparser.add_argument("--errcap", type=int, default=4, help="Width of error bar caps")
    subparser.add_argument("--vertical", action="store_true", help="Stack bars vertically (default True)", default=False)
    if fastq_parser==True:
        subparser.add_argument("--PDB_pt", type=str, default=argparse.SUPPRESS, help="PDB ID (if saved to ~/.config/edms/PDB) or file path for PDB structure file. See edms.dat.pdb.retrieve() or edms uniprot retrieve -h for more information")

    # Figure & layout
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(5,5), help="Figure size formatted as 'width,height'")
    subparser.add_argument("--title", type=str, default="", help="Plot title")
    subparser.add_argument("--title_size", type=int, default=18, help="Font size of the title")
    subparser.add_argument("--title_weight", type=str, default="bold", help="Font weight of the title (e.g., bold, normal)")
    subparser.add_argument("--title_font", type=str, default="Arial", help="Font family for the title")

    # X-axis formatting
    subparser.add_argument("--x_axis", type=str, default="", help="X-axis label")
    subparser.add_argument("--x_axis_size", type=int, default=12, help="Font size for X-axis label")
    subparser.add_argument("--x_axis_weight", type=str, default="bold", help="Font weight for X-axis label")
    subparser.add_argument("--x_axis_font", type=str, default="Arial", help="Font family for X-axis label")
    subparser.add_argument("--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("--x_ticks_size", type=int, default=9, help="Font size for X-axis tick labels")
    subparser.add_argument("--x_ticks_rot", type=int, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")

    # Y-axis formatting
    subparser.add_argument("--y_axis", type=str, default="", help="Y-axis label")
    subparser.add_argument("--y_axis_size", type=int, default=12, help="Font size for Y-axis label")
    subparser.add_argument("--y_axis_weight", type=str, default="bold", help="Font weight for Y-axis label")
    subparser.add_argument("--y_axis_font", type=str, default="Arial", help="Font family for Y-axis label")
    subparser.add_argument("--y_axis_dims", type=parse_tuple_float, default=(0,0), help="Y-axis range as tuple: start,end")
    subparser.add_argument("--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("--y_ticks_size", type=int, default=9, help="Font size for Y-axis tick labels")
    subparser.add_argument("--y_ticks_rot", type=int, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")

    # Legend options
    subparser.add_argument("--legend_title", type=str, default="", help="Legend title text")
    subparser.add_argument("--legend_title_size", type=int, default=12, help="Font size of the legend title")
    subparser.add_argument("--legend_size", type=int, default=12, help="Font size for legend items")
    subparser.add_argument("--legend_bbox_to_anchor", type=parse_tuple_float, default=(1, 1), help="Anchor position for the legend bounding box")
    subparser.add_argument("--legend_loc", type=str, default="upper left", help="Legend location on the plot")
    subparser.add_argument("--legend_ncol", type=int, default=1, help="Number of columns in the legend")
    subparser.add_argument('--legend_columnspacing', type=int, default=argparse.SUPPRESS, help='space between columns in legend; only for html plots')
    subparser.add_argument('--legend_handletextpad', type=float, default=argparse.SUPPRESS, help='space between marker and text in legend; only for html plots')
    subparser.add_argument('--legend_labelspacing', type=float, default=argparse.SUPPRESS, help='vertical space between entries in legend; only for html plots')
    subparser.add_argument('--legend_borderpad', type=float, default=argparse.SUPPRESS, help='padding inside legend box; only for html plots')
    subparser.add_argument('--legend_handlelength', type=float, default=argparse.SUPPRESS, help='marker length in legend; only for html plots')
    subparser.add_argument('--legend_size_html_multiplier', type=float, default=argparse.SUPPRESS, help='legend size multiplier for html plots')

    # Display and formatting
    subparser.add_argument("--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
    subparser.add_argument("--show", action="store_true", help="Show the plot in an interactive window", default=False)
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize and space legend/label values", default=False)

def add_common_plot_vol_args(subparser, fastq_parser=False):
    '''
    add_common_plot_vol_args(subparser): Add common arguments for volcano plot
    '''
    # Required arguments
    subparser.add_argument("--df", type=str, help="Input dataframe file path from edms stat compare", required=True)
    subparser.add_argument("--FC", type=str, help="Fold change column name (X-axis)", required=True)
    subparser.add_argument("--pval", type=str, help="P-value column name (Y-axis)", required=True)

    # Optional data columns
    subparser.add_argument("--stys", type=str, help="Style column name for custom markers")
    subparser.add_argument("--size", type=str, help="Column name used to scale point sizes (Default: -log10('pval'); specify 'false' for no size)")
    subparser.add_argument("--size_dims", type=parse_tuple_float, help="Size range for points formatted as min,max")
    subparser.add_argument("--label", type=str, help="Column containing text labels for points")
    if fastq_parser==False:
        subparser.add_argument("--stys_order", type=str, nargs="+", help="Style column values order")
        subparser.add_argument("--mark_order", type=str, nargs="+", help="Markers order for style column values order")
    
    # Additional annotation data sources
    if fastq_parser==True:
        add_common_fastq_label_args(subparser)
        
    # Thresholds
    subparser.add_argument("--FC_threshold", type=float, default=1, help="Fold change threshold for significance")
    subparser.add_argument("--pval_threshold", type=float, default=1, help="P-value threshold for significance")

    # Output
    subparser.add_argument("--dir", type=str, help="Output directory path", default='./out')
    subparser.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_vol.png')

    # Aesthetics
    subparser.add_argument("--color", type=str, default="lightgray", help="Color for non-significant points")
    subparser.add_argument("--alpha", type=float, default=0.5, help="Transparency for non-significant points")
    subparser.add_argument("--edgecol", type=str, default="black", help="Edge color of points")
    subparser.add_argument("--vertical", action="store_true", help="Use vertical layout for plot", default=False)

    # Figure setup
    subparser.add_argument("--figsize", type=parse_tuple_int, default=(5,5), help="Figure size formatted as 'width,height'")
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
    subparser.add_argument("--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("--x_ticks_size", type=int, default=9, help="Font size for X-axis tick labels")
    subparser.add_argument("--x_ticks_rot", type=int, default=0, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")
    subparser.add_argument("--x_ticks", nargs="+", help="Custom tick values for X-axis")

    # Y-axis settings
    subparser.add_argument("--y_axis", type=str, default="", help="Label for the Y-axis")
    subparser.add_argument("--y_axis_size", type=int, default=12, help="Font size for Y-axis label")
    subparser.add_argument("--y_axis_weight", type=str, default="bold", help="Font weight for Y-axis label")
    subparser.add_argument("--y_axis_font", type=str, default="Arial", help="Font family for Y-axis label")
    subparser.add_argument("--y_axis_dims", type=parse_tuple_float, default=(0, 0), help="Y-axis range formatted as min,max")
    subparser.add_argument("--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("--y_ticks_size", type=int, default=9, help="Font size for Y-axis tick labels")
    subparser.add_argument("--y_ticks_rot", type=int, default=0, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")
    subparser.add_argument("--y_ticks", nargs="+", help="Custom tick values for Y-axis")

    # Legend
    subparser.add_argument("--legend_title", type=str, default="", help="Title for the legend")
    subparser.add_argument("--legend_title_size", type=int, default=12, help="Font size for legend title")
    subparser.add_argument("--legend_size", type=int, default=9, help="Font size for legend items")
    subparser.add_argument("--legend_bbox_to_anchor", type=parse_tuple_float, default=(1, 1), help="Bounding box anchor for legend")
    subparser.add_argument("--legend_loc", type=str, default="upper left", help="Legend location on the plot")
    subparser.add_argument("--legend_ncol", type=int, default=1, help="Number of columns in the legend")
    subparser.add_argument('--legend_columnspacing', type=int, default=argparse.SUPPRESS, help='space between columns in legend; only for html plots')
    subparser.add_argument('--legend_handletextpad', type=float, default=argparse.SUPPRESS, help='space between marker and text in legend; only for html plots')
    subparser.add_argument('--legend_labelspacing', type=float, default=argparse.SUPPRESS, help='vertical space between entries in legend; only for html plots')
    subparser.add_argument('--legend_borderpad', type=float, default=argparse.SUPPRESS, help='padding inside legend box; only for html plots')
    subparser.add_argument('--legend_handlelength', type=float, default=argparse.SUPPRESS, help='marker length in legend; only for html plots')
    subparser.add_argument('--legend_size_html_multiplier', type=float, default=argparse.SUPPRESS, help='legend size multiplier for html plots')

    # Boolean switches
    subparser.add_argument("--dont_display_legend", action="store_false", help="Don't display legend on plot", default=True)
    subparser.add_argument("--display_labels", type=str, nargs="+", help="Display labels for values if label column specified (Options: 'FC & p-value', 'FC', 'p-value', 'NS', 'all', or ['label1', 'label2', ..., 'labeln'])", default=["FC & p-value"])
    subparser.add_argument("--dont_display_axis", dest='display_axis', action="store_false", default=True, help="Display x- and y-axis lines (Default: True)")
    subparser.add_argument("--display_lines", action="store_true", help="Display lines for threshold (Default: False)", default=False)
    subparser.add_argument("--return_df", action="store_true", help="Return annotated DataFrame after plotting", default=False)
    subparser.add_argument("--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
    subparser.add_argument("--show", action="store_true", help="Show the plot in an interactive window", default=False)
    subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize and space labels/legend items", default=False)

def add_subparser(subparsers, formatter_class=None):
    """
    add_subparser(): Attach all gen-related subparsers to the top-level CLI.

    Parameters:
    subparsers (argparse._SubParsersAction): The subparsers object to attach the gen subparsers to.
    formatter_class (type, optional): The formatter class to use for the subparsers.
    """
    '''
    edms.gen.plot:
    - scat(): creates scatter plot related graphs
    - cat(): creates categorical graphs
    - dist(): creates distribution graphs
    - heat(): creates heatmap graphs
    - stack(): creates stacked bar plot
    - vol(): creates volcano plot
    '''
    parser_plot = subparsers.add_parser("plot", help="Generate scatter, category, distribution, heatmap, stacked bar, and volcano plots", description="Generate scatter, category, distribution, heatmap, stacked bar, and volcano plots", formatter_class=formatter_class)
    subparsers_plot = parser_plot.add_subparsers(dest="typ")

    # scat(): Creates scatter plot related graphs (scat, line, line_scat)
    parser_plot_type_scat = subparsers_plot.add_parser("scat", help="Create scatter plot", description="Create scatter plot", formatter_class=formatter_class)
    parser_plot_type_line = subparsers_plot.add_parser("line", help="Create line plot", description="Create line plot", formatter_class=formatter_class)
    parser_plot_type_line_scat = subparsers_plot.add_parser("line_scat", help="Create scatter + line plot", description="Create scatter + line plot", formatter_class=formatter_class)

    for parser_plot_scat in [parser_plot_type_scat, parser_plot_type_line, parser_plot_type_line_scat]:
        add_common_plot_scat_args(parser_plot_scat)
        parser_plot_scat.set_defaults(func=p.scat)

    # cat(): Creates categorical graphs (bar, box, violin, swarm, strip, point, count, bar_swarm, box_swarm, violin_swarm)
    parser_plot_type_bar = subparsers_plot.add_parser("bar", help="Create bar plot", description="Create bar plot", formatter_class=formatter_class)
    parser_plot_type_box = subparsers_plot.add_parser("box", help="Create box plot", description="Create box plot", formatter_class=formatter_class)
    parser_plot_type_violin = subparsers_plot.add_parser("violin", help="Create violin plot", description="Create violin plot", formatter_class=formatter_class)
    parser_plot_type_swarm = subparsers_plot.add_parser("swarm", help="Create swarm plot", description="Create swarm plot", formatter_class=formatter_class)
    parser_plot_type_strip = subparsers_plot.add_parser("strip", help="Create strip plot", description="Create strip plot", formatter_class=formatter_class)
    parser_plot_type_point = subparsers_plot.add_parser("point", help="Create point plot", description="Create point plot", formatter_class=formatter_class)
    parser_plot_type_count = subparsers_plot.add_parser("count", help="Create count plot", description="Create count plot", formatter_class=formatter_class)
    parser_plot_type_bar_swarm = subparsers_plot.add_parser("bar_swarm", help="Create bar + swarm plot", description="Create bar + swarm plot", formatter_class=formatter_class)
    parser_plot_type_box_swarm = subparsers_plot.add_parser("box_swarm", help="Create box + swarm plot", description="Create box + swarm plot", formatter_class=formatter_class)
    parser_plot_type_violin_swarm = subparsers_plot.add_parser("violin_swarm", help="Create violin + swarm plot", description="Create violin + swarm plot", formatter_class=formatter_class)

    for parser_plot_cat in [parser_plot_type_bar, parser_plot_type_box, parser_plot_type_violin, parser_plot_type_swarm, parser_plot_type_strip, parser_plot_type_point, parser_plot_type_count, parser_plot_type_bar_swarm, parser_plot_type_box_swarm, parser_plot_type_violin_swarm]:
        add_common_plot_cat_args(parser_plot_cat)
        parser_plot_cat.set_defaults(func=p.cat)

    # dist(): Creates distribution graphs (hist, kde, hist_kde, rid)
    parser_plot_type_hist = subparsers_plot.add_parser("hist", help="Create histogram plot", description="Create histogram plot", formatter_class=formatter_class)
    parser_plot_type_kde = subparsers_plot.add_parser("kde", help="Create density plot", description="Create density plot", formatter_class=formatter_class)
    parser_plot_type_hist_kde = subparsers_plot.add_parser("hist_kde", help="Create histogram + density plot", description="Create histogram + density plot", formatter_class=formatter_class)
    parser_plot_type_rid = subparsers_plot.add_parser("rid", help="Create ridge plot", description="Create ridge plot", formatter_class=formatter_class)

    for parser_plot_dist in [parser_plot_type_hist, parser_plot_type_kde, parser_plot_type_hist_kde, parser_plot_type_rid]:
        add_common_plot_dist_args(parser_plot_dist)
        parser_plot_dist.set_defaults(func=p.dist)

    # heat(): Creates heatmap graphs
    parser_plot_type_heat = subparsers_plot.add_parser("heat", help="Create heatmap plot", description="Create heatmap plot", formatter_class=formatter_class)
    add_common_plot_heat_args(parser_plot_type_heat)
    parser_plot_type_heat.set_defaults(func=p.heat)
    
    # stack(): Creates stacked bar plot
    parser_plot_type_stack = subparsers_plot.add_parser("stack", help="Create stacked bar plot", description="Create stacked bar plot", formatter_class=formatter_class)
    add_common_plot_stack_args(parser_plot_type_stack)
    parser_plot_type_stack.set_defaults(func=p.stack)

    # vol(): Creates volcano plot
    parser_plot_type_vol = subparsers_plot.add_parser("vol", help="Create volcano plot", description="Create volcano plot", formatter_class=formatter_class)
    add_common_plot_vol_args(parser_plot_type_vol)
    parser_plot_type_vol.set_defaults(func=p.vol)

    '''
    edms.gen.stat:
    - describe(): returns descriptive statistics for numerical columns in a DataFrame
    - difference(): computes the appropriate statistical test(s) and returns the p-value(s)
    - correlation(): returns a correlation matrix
    - compare(): computes FC, pval, and log transformations relative to a specified condition
    '''
    parser_stat = subparsers.add_parser("stat", help="Statistics", description="Statistics", formatter_class=formatter_class)
    subparsers_stat = parser_stat.add_subparsers()
    
    # describe(): returns descriptive statistics for numerical columns in a DataFrame
    parser_stat_describe = subparsers_stat.add_parser("describe", help="Compute descriptive statistics", description="Compute descriptive statistics", formatter_class=formatter_class)

    parser_stat_describe.add_argument("--df", type=str, help="Input file path", required=True)

    parser_stat_describe.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_stat_describe.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_descriptive.csv')
    
    parser_stat_describe.add_argument("--cols", nargs="+", help="List of numerical columns to describe")
    parser_stat_describe.add_argument("--group", type=str, help="Column name to group by")
    
    parser_stat_describe.set_defaults(func=st.describe)

    # difference(): computes the appropriate statistical test(s) and returns the p-value(s)
    parser_stat_difference = subparsers_stat.add_parser("difference", help="Compute statistical difference between groups", description="Compute statistical difference between groups", formatter_class=formatter_class)

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
    parser_stat_correlation = subparsers_stat.add_parser("correlation", help="Compute correlation matrix", description="Compute correlation matrix", formatter_class=formatter_class)

    parser_stat_correlation.add_argument("--df", type=str, help="Input file path",required=True)
    parser_stat_correlation.add_argument("--var_cols", nargs="+", help="List of 2 variable columns for tidy format")
    parser_stat_correlation.add_argument("--value_cols", nargs="+", help="List of numerical columns to correlate")
    parser_stat_correlation.add_argument("--method", type=str, default="pearson", choices=["pearson", "spearman", "kendall"],
                                         help="Correlation method to use (Default: pearson)")
    parser_stat_correlation.add_argument("--numeric_only", action="store_true", help="Only use numeric columns (Default: True)")
    parser_stat_correlation.add_argument("--no_plot", dest="plot", action="store_false", help="Don't generate correlation matrix plot", default=True)
    parser_stat_correlation.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_stat_correlation.add_argument("--file_data", type=str, help="Output data file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_correlation.csv')
    parser_stat_correlation.add_argument("--file_plot", type=str, help="Output plot file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_correlation.pdf')
    add_common_plot_heat_args(parser_stat_correlation, stat_parser=True)

    parser_stat_correlation.set_defaults(func=st.correlation)

    # compare(): computes FC, pval, and log transformations relative to a specified condition
    parser_stat_compare = subparsers_stat.add_parser("compare", help="Compare conditions using FC, p-values, and log transforms", description="Compare conditions using FC, p-values, and log transforms", formatter_class=formatter_class)

    parser_stat_compare.add_argument("--df", type=str, help="Input file path",required=True)
    parser_stat_compare.add_argument("--sample", type=str, help="Sample column name",required=True)
    parser_stat_compare.add_argument("--cond", type=str, help="Condition column name",required=True)
    parser_stat_compare.add_argument("--cond_comp", type=str, help="Condition to compare against",required=True)
    parser_stat_compare.add_argument("--var", type=str, help="Variable column name",required=True)
    parser_stat_compare.add_argument("--count", type=str, help="Count column name",required=True)

    parser_stat_compare.add_argument("--pseudocount", type=int, default=1, help="Pseudocount to avoid log(0) or divide-by-zero errors")
    parser_stat_compare.add_argument("--alternative", type=str, default="two-sided", choices=["two-sided", "less", "greater"], help="Alternative hypothesis for Fisher's exact test (Default: two-sided)")
    parser_stat_compare.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_stat_compare.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_compare.csv')
    parser_stat_compare.add_argument("--verbose", action="store_true", help="Print progress to console", default=False)

    parser_stat_compare.set_defaults(func=st.compare)

    # odds_ratio(): computes odds ratio relative to a specified condition (OR = (A/B)/(C/D))
    parser_stat_odds_ratio = subparsers_stat.add_parser("odds_ratio", help="Computes odds ratios relative to a specified condition & variable (e.g., unedited & WT)", description="Compute odds ratio between conditions", formatter_class=formatter_class)

    parser_stat_odds_ratio.add_argument("--df", type=str, help="Input file path",required=True)
    parser_stat_odds_ratio.add_argument("--cond", type=str, help="Condition column name",required=True)
    parser_stat_odds_ratio.add_argument("--cond_comp", type=str, help="Condition for comparison group",required=True)
    parser_stat_odds_ratio.add_argument("--var", type=str, help="Variable column name",required=True)
    parser_stat_odds_ratio.add_argument("--var_comp", type=str, help="Variable name for comparison (e.g., WT)",required=True)
    parser_stat_odds_ratio.add_argument("--count", type=str, help="Count column name",required=True)
    
    parser_stat_odds_ratio.add_argument("--pseudocount", type=int, default=1, help="Pseudocount to avoid /0 (Default: 1)")
    parser_stat_odds_ratio.add_argument("--alternative", type=str, default="two-sided", choices=["two-sided", "less", "greater"], help="Alternative hypothesis for Fisher's exact test (Default: two-sided)")
    parser_stat_odds_ratio.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_stat_odds_ratio.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_odds_ratio.csv')
    parser_stat_odds_ratio.add_argument("--verbose", action="store_true", help="Print progress to console", default=False)

    parser_stat_odds_ratio.set_defaults(func=st.odds_ratio)

    '''
    edms.gen.io:
    - in_subs(): moves all files with a given suffix into subfolders named after the files (excluding the suffix).
    - out_subs(): recursively moves all files from subdirectories into the parent directory and delete the emptied subdirectories.
    - create_sh(): creates a shell script with SLURM job submission parameters for Harvard FASRC cluster.
    - split_R1_R2(): split paired reads into new R1 and R2 subdirectories at the parent directory
    - excel_csvs(): exports excel file to .csv files in specified directory  
    '''
    parser_io = subparsers.add_parser("io", help="Input/Output", formatter_class=formatter_class)
    subparsers_io = parser_io.add_subparsers()
    
    # Create subparsers for commands
    parser_io_in_subs = subparsers_io.add_parser("in_subs", help="*No FASRC* Moves all files with a given suffix into subfolders named after the files (excluding the suffix)", description="*No FASRC* Moves all files with a given suffix into subfolders named after the files (excluding the suffix)", formatter_class=formatter_class)
    parser_io_out_subs = subparsers_io.add_parser("out_subs", help="*No FASRC* Delete subdirectories and move their files to the parent directory", description="*No FASRC* Delete subdirectories and move their files to the parent directory", formatter_class=formatter_class)
    parser_io_create_sh = subparsers_io.add_parser("create_sh", help='Generate SLURM shell script for Harvard FASRC cluster.', description='Generate SLURM shell script for Harvard FASRC cluster.', formatter_class=formatter_class)
    parser_io_split_R1_R2 = subparsers_io.add_parser("split_R1_R2", help='*No FASRC* Split paired reads into new R1 and R2 subdirectories at the parent directory.', description='*No FASRC* Split paired reads into new R1 and R2 subdirectories at the parent directory.', formatter_class=formatter_class)
    parser_io_excel_csvs = subparsers_io.add_parser("excel_csvs", help='Exports excel file to .csv files in specified directory.', description='Exports excel file to .csv files in specified directory.', formatter_class=formatter_class)

    # Add common arguments
    for parser_io_common in [parser_io_in_subs,parser_io_out_subs,parser_io_split_R1_R2]:
        parser_io_common.add_argument("--dir", help="Path to parent directory", type=str, default='.')
    
    # in_subs() arguments
    parser_io_in_subs.add_argument("--suf", help="File suffix (e.g., '.txt', '.csv') to filter files.", type=str, required=True) 
    
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
    edms.gen.com:
    - access(): make all files and subdirectories accessible on Harvard FASRC
    - smaller_fastq(): create new subdirectory containing fastqs with the # of reads limited
    - create_export_var(): create a persistent environment variable by adding it to the user's shell config.
    - view_export_vars(): View the current export variables in the user's shell config.
    '''
    parser_com = subparsers.add_parser("com", help="Command Line Interaction", description="Command Line Interaction", formatter_class=formatter_class)
    subparsers_com = parser_com.add_subparsers()
    
    # Create subparsers for commands
    parser_com_access = subparsers_com.add_parser("access", help="Make all files and subdirectories accessible on Harvard FASRC", description="Make all files and subdirectories accessible on Harvard FASRC", formatter_class=formatter_class)
    parser_com_smaller_fastq = subparsers_com.add_parser("smaller_fastq", help="Create new subdirectory containing fastqs with the # of reads limited", description="Create new subdirectory containing fastqs with the # of reads limited", formatter_class=formatter_class)
    parser_com_create_export_var = subparsers_com.add_parser("create_export_var", help="Create a persistent export variable by adding it to the user's shell config.", description="Create a persistent export variable by adding it to the user's shell config.", formatter_class=formatter_class)
    parser_com_view_export_vars = subparsers_com.add_parser("view_export_vars", help="View the current export variables in the user's shell config.", description="View the current export variables in the user's shell config.", formatter_class=formatter_class)

    # Add common arguments
    for parser_com_common in [parser_com_access, parser_com_smaller_fastq]:
        parser_com_common.add_argument("--pt", help="Path to parent directory", type=str, default='.')
    
    # Smaller_fastq arguments
    parser_com_smaller_fastq.add_argument("--reads", help="# of reads per fastq file", type=int, default='100000') 
    parser_com_smaller_fastq.add_argument("--suf", help="Fastq file suffix", type=int, default=".fastq.gz") 
    
    # create_export_var arguments
    parser_com_create_export_var.add_argument("--name", help="Name of the environment variable (e.g., MYPROJ)", required=True)
    parser_com_create_export_var.add_argument("--pt", help="Path the variable should point to (e.g., ~/projects/myproj)", required=True)
    parser_com_create_export_var.add_argument("--shell", choices=["bash", "zsh"], default=argparse.SUPPRESS, help="Shell type)")
    
    # view_export_var arguments
    parser_com_view_export_vars.add_argument("--shell", choices=["bash", "zsh"], default=argparse.SUPPRESS, help="Shell type")

    # set default functions
    parser_com_access.set_defaults(func=com.access)
    parser_com_smaller_fastq.set_defaults(func=com.smaller_fastq)
    parser_com_create_export_var.set_defaults(func=com.create_export_var)
    parser_com_view_export_vars.set_defaults(func=com.view_export_vars)

    '''
    edms.gen.html:
    - make_html_index(): Create an index HTML that links to other HTML files in `dir`.
        - Uses <title> from each HTML file if available; falls back to stem/filename.
        - Makes titles into buttons.
        - Optionally embeds a preview iframe that updates when you click a button.
        - Displays plot links in a responsive grid; `grid_cols` controls the default column count.
    '''
    parser_html = subparsers.add_parser("html", help="HTML Index Creation", description="HTML Index Creation", formatter_class=formatter_class)
    
    parser_html.add_argument("--dir", type=str, help="Directory containing HTML files to index", default=".")
    parser_html.add_argument("--file", type=str, help="Output HTML index file name", default="index.html")
    parser_html.add_argument("--recursive", action="store_true", help="Recursively search subdirectories for HTML files", default=False)
    parser_html.add_argument("--exclude", type=str, nargs="+", help="List of filenames to exclude (case insensitive)", default=[])
    parser_html.add_argument("--sort", type=str, choices=["title", "name", "mtime"], help="Sort HTML files by 'title', 'name', or 'mtime' (modification time)", default="title")
    parser_html.add_argument("--no_preview", dest="preview", action="store_false", help="Don't include an iframe preview panel in the index", default=True)
    parser_html.add_argument("--grid_cols", type=int, help="Number of columns in the responsive grid layout", default=3)
    parser_html.add_argument("--image_types", type=str, nargs="+", help="List of image file extensions to include (e.g. .png .jpg .gif). If not specified, only .html files are included.", default=None)
    parser_html.add_argument("--preview_height_px", type=int, help="Height of the preview iframe in pixels", default=900)
    parser_html.add_argument("--icon", type=str, help="Name of the SVG icon file (without .svg) to use as favicon", default="python")
    
    parser_html.set_defaults(func=ht.make_html_index)