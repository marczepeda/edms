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
    subparser.add_argument("-L", "--label_size", type=int, help="Label font size (Default: 16)", default=16)
    subparser.add_argument("-N", "--no_label_info", dest='label_info', action="store_false", help="Don't include additional info for labels if .html plot (Default: True)", default=True)
    subparser.add_argument("-aa", "--aa_properties", nargs="+", help="Use aa_properties to format labels (Options: hydrophobicity, polarity, charge, vdw_volume, pKa_C_term, pKa_N_term, pKa_side_chain)", default=argparse.SUPPRESS)
    subparser.add_argument("-cBP", "--cBioPortal", type=str, help="Gene name (if saved to ~/.config/edms/cBioPortal_mutations) or file path for cBioPortal mutation data processed through edms.dat.cBioPortal.mutations()", default=argparse.SUPPRESS)
    subparser.add_argument("-UP", "--UniProt", type=str, help="UniProt accession (if saved to ~/.config/edms/UniProt) or file path for UniProt flat file. See edms.dat.uniprot.retrieve_flat_file() or edms uniprot retrieve -h for more information.", default=argparse.SUPPRESS)
    subparser.add_argument("-PSP", "--PhosphoSitePlus", type=str, help="UniProt accession", default=argparse.SUPPRESS)
    subparser.add_argument("-PDBc", "--PDB_contacts", type=str, help="PDB ID (if saved to ~/.config/edms/PDB) or file path for PDB structure file. See edms.dat.pdb.retrieve() or edms uniprot retrieve -h for more information.", default=argparse.SUPPRESS)
    subparser.add_argument("-PDBn", "--PDB_neighbors", type=str, help="PDB ID (if saved to ~/.config/edms/PDB) or file path for PDB structure file. See edms.dat.pdb.retrieve() or edms uniprot retrieve -h for more information.", default=argparse.SUPPRESS)
# Plot subparser methods
def add_common_plot_scat_args(subparser, fastq_torn_parser=False, fastq_corr_parser=False):
    '''
    add_common_plot_scat_args(subparser): Add common arguments for scatter plot related graphs
    '''
    # scat(): Required arguments
    subparser.add_argument("-i", "--df", help="Input dataframe file path", type=str, required=True)
    if fastq_torn_parser == False and fastq_corr_parser == False:
        subparser.add_argument("-x", "--x", help="X-axis column", type=str, required=True)
        subparser.add_argument("-y", "--y", help="Y-axis column", type=str, required=True)
    else:
        if fastq_corr_parser == True:
            subparser.add_argument("-cc", "--cond_col", help="condition column name for comparison", type=str, required=True)
            subparser.add_argument("-cv", "--cond_vals", nargs="+", help="two condition values for comparison (x and y-axis)", type=str, required=True)
        subparser.add_argument("-FC", "--FC", help="Fold change column name (Y-axis)", type=str, required=True)
        subparser.add_argument("-pval", "--pval", help="p-value column name (size column if not specified)", type=str, required=True)
        
    # Optional core arguments
    if fastq_torn_parser == False and fastq_corr_parser == False:
        subparser.add_argument("-c", "--cols", type=str, help="Color column name")
        subparser.add_argument("-co", "--cols_ord", nargs="+", help="Column order (list of values)")
        subparser.add_argument("-ce", "--cols_exclude", nargs="+", help="Columns to exclude from coloring")
        subparser.add_argument("-st", "--stys", type=str, help="Style column name")
        subparser.add_argument("-so", "--stys_order", nargs="+", help="Style order (list of values)")
        subparser.add_argument("-mo", "--mark_order", nargs="+", help="Marker order (list of marker styles)")
    else:
        subparser.add_argument("-si", "--size", type=str, help="Column name used to scale point sizes (Default: -log10('pval'); specify 'false' for no size)")
        subparser.add_argument("-sd", "--size_dims", type=parse_tuple_float, help="Size range for points formatted as min,max")
        subparser.add_argument("-zc", "--z_col", type=str, help="Column name for Z-score normalization (Default: None)", default=argparse.SUPPRESS)
        subparser.add_argument("-zv", "--z_var", type=str, help="Column value for Z-score normalization (Default: None)", default=argparse.SUPPRESS)
        if fastq_corr_parser == True:
            subparser.add_argument("-m", "--method", help="Correlation method (Default: 'pearson')", choices=['pearson', 'spearman', 'kendall'], type=str, default='pearson')
            subparser.add_argument("-nw", "--not_weighted", dest='weighted', action='store_false', help="Weighted correlation by size column (Default: True)", default=True)

    subparser.add_argument("-l", "--label", type=str, help="Column name for point labels; static text for images, interactive tooltips for HTML")

    # Additional annotation data sources
    if fastq_torn_parser==True or fastq_corr_parser==True:
        add_common_fastq_label_args(subparser)
        if fastq_torn_parser==True:
            subparser.add_argument("-ss_h", "--ss_h", type=int, help="Height for secondary structure in the plot (Default: autogenerate)", default=argparse.SUPPRESS)
            subparser.add_argument("-ss_y", "--ss_y", type=int, help="Y position for secondary structure track in the plot (Default: autogenerate)", default=argparse.SUPPRESS)

    subparser.add_argument("-o", "--dir", help="Output directory path", type=str, default='./out')
    subparser.add_argument("-f", "--file", help="Output file name", type=str, required=False, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_scat.png')
    if fastq_torn_parser == False and fastq_corr_parser == False:
        subparser.add_argument("-pc", "--palette_or_cmap", type=str, default="colorblind", help="Seaborn palette or matplotlib colormap")
        subparser.add_argument('-a', "--alpha", type=float, default=1, help="Alpha (transparency) for scatter points (0 to 1)")
    subparser.add_argument("-ec", "--edgecol", type=str, default="black", help="Edge color for scatter points")

    # Figure appearance
    subparser.add_argument("-fs", "--figsize", type=parse_tuple_int, default=(5,5), help="Figure size as a tuple: width,height")
    subparser.add_argument("-t", "--title", type=str, default="", help="Plot title")
    subparser.add_argument("-ts", "--title_size", type=int, default=18, help="Plot title font size")
    subparser.add_argument("-tw", "--title_weight", type=str, default="bold", help="Plot title font weight (e.g., bold, normal)")
    subparser.add_argument("-tf", "--title_font", type=str, default="Arial", help="Font family for the title")
    # X-axis settings
    subparser.add_argument("-xa", "--x_axis", type=str, default="", help="X-axis label")
    subparser.add_argument("-xas", "--x_axis_size", type=int, default=12, help="X-axis label font size")
    subparser.add_argument("-xaw", "--x_axis_weight", type=str, default="bold", help="X-axis label font weight (e.g., bold)")
    subparser.add_argument("-xaf", "--x_axis_font", type=str, default="Arial", help="X-axis label font family")
    subparser.add_argument("-xasc", "--x_axis_scale", type=str, default="linear", help="X-axis scale: linear, log, etc.")
    subparser.add_argument("-xad", "--x_axis_dims", type=parse_tuple_int, default=(0,0), help="X-axis range as a tuple: start,end")
    subparser.add_argument("-xap", "--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("-xts", "--x_ticks_size", type=int, default=9, help="X-axis tick labels font size")
    subparser.add_argument("-xtr", "--x_ticks_rot", type=int, default=0, help="Rotation angle of X-axis tick labels")
    subparser.add_argument("-xtf", "--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")
    subparser.add_argument("-xt", "--x_ticks", nargs="+", help="Specific tick values for X-axis")

    # Y-axis settings
    subparser.add_argument("-ya", "--y_axis", type=str, default="", help="Y-axis label")
    subparser.add_argument("-yas", "--y_axis_size", type=int, default=12, help="Y-axis label font size")
    subparser.add_argument("-yaw", "--y_axis_weight", type=str, default="bold", help="Y-axis label font weight")
    subparser.add_argument("-yaf", "--y_axis_font", type=str, default="Arial", help="Y-axis label font family")
    subparser.add_argument("-yasc", "--y_axis_scale", type=str, default="linear", help="Y-axis scale: linear, log, etc.")
    subparser.add_argument("-yad", "--y_axis_dims", type=parse_tuple_int, default=(0,0), help="Y-axis range as a tuple: start,end")
    subparser.add_argument("-yap", "--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("-yts", "--y_ticks_size", type=int, default=9, help="Y-axis tick labels font size")
    subparser.add_argument("-ytr", "--y_ticks_rot", type=int, default=0, help="Rotation angle of Y-axis tick labels")
    subparser.add_argument("-ytf", "--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")
    subparser.add_argument("-yt", "--y_ticks", nargs="+", help="Specific tick values for Y-axis")

    # Legend settings
    subparser.add_argument("-lt", "--legend_title", type=str, default="", help="Legend title")
    subparser.add_argument("-lts", "--legend_title_size", type=int, default=12, help="Legend title font size")
    subparser.add_argument("-ls", "--legend_size", type=int, default=9, help="Legend font size")
    subparser.add_argument("-lba", "--legend_bbox_to_anchor", type=parse_tuple_float, default=(1,1), help="Bounding box anchor position for legend")
    subparser.add_argument("-ll", "--legend_loc", type=str, default="upper left", help="Location of the legend in the plot")
    subparser.add_argument("-li", "--legend_items", type=parse_tuple_int, default=(0,0), help="Legend item count as a tuple (used for layout)")
    subparser.add_argument("-ln", "--legend_ncol", type=int, default=1, help="Number of columns in legend")
    subparser.add_argument('-lcs', "--legend_columnspacing", type=int, default=argparse.SUPPRESS, help='space between columns in legend; only for html plots')
    subparser.add_argument('-lhtp', "--legend_handletextpad", type=float, default=argparse.SUPPRESS, help='space between marker and text in legend; only for html plots')
    subparser.add_argument('-lls', "--legend_labelspacing", type=float, default=argparse.SUPPRESS, help='vertical space between entries in legend; only for html plots')
    subparser.add_argument('-lbp', "--legend_borderpad", type=float, default=argparse.SUPPRESS, help='padding inside legend box; only for html plots')
    subparser.add_argument('-lhl', "--legend_handlelength", type=float, default=argparse.SUPPRESS, help='marker length in legend; only for html plots')
    subparser.add_argument('-lshm', "--legend_size_html_multiplier", type=float, default=argparse.SUPPRESS, help='legend size multiplier for html plots')

    # Display and formatting
    subparser.add_argument("-d", "--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
    subparser.add_argument("-s", "--show", action="store_true", help="Show the plot", default=False)
    subparser.add_argument("-sc", "--space_capitalize", action="store_true", help="Capitalize label/legend strings and replace underscores with spaces")
    if fastq_torn_parser == True or fastq_corr_parser == True:
        subparser.add_argument("-ddlg", "--dont_display_legend", dest='display_legend', action="store_false", default=True, help="Display legend on plot (Default: True)")
        subparser.add_argument("-ddla", "--dont_display_labels", dest='display_labels', action="store_false", default=True, help="Display labels for significant values (Default: True)")
        subparser.add_argument("-dda", "--dont_display_axis", dest='display_axis', action="store_false", default=True, help="Display x- and y-axis lines (Default: True)")

def add_common_plot_cat_args(subparser, fastq_parser=False):
    '''
    add_common_plot_cat_args(subparser): Add common arguments for categorical graphs
    '''
    # cat(): Required arguments
    if fastq_parser == True:
        subparser.add_argument("-type", "--type", dest="typ", help="Type of category plot", type=str, required=True, choices=['bar', 'box', 'violin', 'strip', 'swarm', 'point', 'count', 'bar_strip', 'box_strip', 'violin_strip', 'bar_swarm', 'box_swarm', 'violin_swarm'])
    subparser.add_argument("-i", "--df", help="Input dataframe file path", type=str, required=True)

    # Optional core arguments
    subparser.add_argument("-x", "--x", help="X-axis column name", type=str, default="")
    subparser.add_argument("-y", "--y", help="Y-axis column name", type=str, default="")
    subparser.add_argument("-co", "--cats_ord", nargs="+", help="Category column values order (x- or y-axis)")
    subparser.add_argument("-ce", "--cats_exclude", nargs="+", help="Category column values exclude (x- or y-axis)")
    subparser.add_argument("-cl", "--cols", type=str, help="Color column name for grouping")
    subparser.add_argument("-clo", "--cols_ord", nargs="+", help="Color column values order")
    subparser.add_argument("-cle", "--cols_exclude", nargs="+", help="Color column values to exclude")
    if fastq_parser == True:
        subparser.add_argument("-PDB_pt", "--PDB_pt", type=str, help="PDB ID (if saved to ~/.config/edms/PDB) or file path for PDB structure file. See edms.dat.pdb.retrieve() or edms uniprot retrieve -h for more information", default=argparse.SUPPRESS)

    subparser.add_argument("-o", "--dir", type=str, help="Output directory", default='./out')
    subparser.add_argument("-f", "--file", type=str, help="Output filename", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_cat.png')
    subparser.add_argument("-pc", "--palette_or_cmap", type=str, default="colorblind", help="Seaborn color palette or matplotlib colormap")
    subparser.add_argument('-a', "--alpha", type=float, default=1, help="Alpha (transparency) for scatter points (0 to 1)")
    subparser.add_argument("-do", "--dodge", action='store_false', default=False, help="Separate points by color category")
    subparser.add_argument("-nj", "--no_jitter", dest="jitter", action='store_false', default=True, help="Don't add jitter for points in strip plots")
    subparser.add_argument("-ms", "--size", type=int, default=5, help="Marker size for scatter points")
    subparser.add_argument("-ec", "--edgecol", type=str, default="black", help="Edge color for markers")

    # Error bar and style options
    subparser.add_argument("-lw", "--lw", type=int, default=1, help="Line width for plot edges")
    subparser.add_argument("-eb", "--errorbar", type=str, default="sd", help="Error bar type: sd (standard deviation), etc.")
    subparser.add_argument("-ew", "--errwid", type=float, default=1, help="Width of the error bars")
    subparser.add_argument("-ecap", "--errcap", type=float, default=0.1, help="Cap size on error bars")

    # Figure appearance
    subparser.add_argument("-fs", "--figsize", type=parse_tuple_int, default=(5,5), help="Figure size formatted 'width,height'")
    subparser.add_argument("-t", "--title", type=str, default="", help="Plot title text")
    subparser.add_argument("-ts", "--title_size", type=int, default=18, help="Font size of the plot title")
    subparser.add_argument("-tw", "--title_weight", type=str, default="bold", help="Font weight of the plot title (e.g., bold, normal)")
    subparser.add_argument("-tf", "--title_font", type=str, default="Arial", help="Font family for the plot title")
    # X-axis settings
    subparser.add_argument("-xa", "--x_axis", type=str, default="", help="X-axis label text")
    subparser.add_argument("-xas", "--x_axis_size", type=int, default=12, help="Font size for the X-axis label")
    subparser.add_argument("-xaw", "--x_axis_weight", type=str, default="bold", help="Font weight for the X-axis label")
    subparser.add_argument("-xaf", "--x_axis_font", type=str, default="Arial", help="Font family for the X-axis label")
    subparser.add_argument("-xasc", "--x_axis_scale", type=str, default="linear", help="Scale of X-axis (e.g., linear, log)")
    subparser.add_argument("-xad", "--x_axis_dims", type=parse_tuple_float, default=(0, 0), help="X-axis range as tuple: start,end")
    subparser.add_argument("-xap", "--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("-xts", "--x_ticks_size", type=int, default=9, help="Font size for X-axis tick labels")
    subparser.add_argument("-xtr", "--x_ticks_rot", type=int, default=0, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("-xtf", "--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")
    subparser.add_argument("-xt", "--x_ticks", nargs="+", help="Explicit tick values for X-axis")

    # Y-axis settings
    subparser.add_argument("-ya", "--y_axis", type=str, default="", help="Y-axis label text")
    subparser.add_argument("-yas", "--y_axis_size", type=int, default=12, help="Font size for the Y-axis label")
    subparser.add_argument("-yaw", "--y_axis_weight", type=str, default="bold", help="Font weight for the Y-axis label")
    subparser.add_argument("-yaf", "--y_axis_font", type=str, default="Arial", help="Font family for the Y-axis label")
    subparser.add_argument("-yasc", "--y_axis_scale", type=str, default="linear", help="Scale of Y-axis (e.g., linear, log)")
    subparser.add_argument("-yad", "--y_axis_dims", type=parse_tuple_float, default=(0, 0), help="Y-axis range as tuple: start,end")
    subparser.add_argument("-yap", "--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("-yts", "--y_ticks_size", type=int, default=9, help="Font size for Y-axis tick labels")
    subparser.add_argument("-ytr", "--y_ticks_rot", type=int, default=0, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("-ytf", "--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")
    subparser.add_argument("-yt", "--y_ticks", nargs="+", help="Explicit tick values for Y-axis")

    # Legend settings
    subparser.add_argument("-lt", "--legend_title", type=str, default="", help="Title for the legend")
    subparser.add_argument("-lts", "--legend_title_size", type=int, default=12, help="Font size for the legend title")
    subparser.add_argument("-ls", "--legend_size", type=int, default=9, help="Font size for legend items")
    subparser.add_argument("-lba", "--legend_bbox_to_anchor", type=parse_tuple_float, default=(1, 1), help="Anchor position of the legend bounding box")
    subparser.add_argument("-ll", "--legend_loc", type=str, default="upper left", help="Location of the legend on the plot")
    subparser.add_argument("-li", "--legend_items", type=parse_tuple_int, default=(0, 0), help="Tuple for legend item layout")
    subparser.add_argument("-ln", "--legend_ncol", type=int, default=1, help="Number of columns in the legend")
    subparser.add_argument('-lcs', "--legend_columnspacing", type=int, default=argparse.SUPPRESS, help='space between columns in legend; only for html plots')
    subparser.add_argument('-lhp', "--legend_handletextpad", type=float, default=argparse.SUPPRESS, help='space between marker and text in legend; only for html plots')
    subparser.add_argument('-lls', "--legend_labelspacing", type=float, default=argparse.SUPPRESS, help='vertical space between entries in legend; only for html plots')
    subparser.add_argument('-lbp', "--legend_borderpad", type=float, default=argparse.SUPPRESS, help='padding inside legend box; only for html plots')
    subparser.add_argument('-lhl', "--legend_handlelength", type=float, default=argparse.SUPPRESS, help='marker length in legend; only for html plots')
    subparser.add_argument('-lshm', "--legend_size_html_multiplier", type=float, default=argparse.SUPPRESS, help='legend size multiplier for html plots')

    # Display and formatting
    subparser.add_argument("-d", "--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
    subparser.add_argument("-s", "--show", action="store_true", help="Show the plot in a window", default=False)
    subparser.add_argument("-sc", "--space_capitalize", action="store_true", help="Capitalize labels and replace underscores with spaces")

def add_common_plot_dist_args(subparser):
    '''
    add_common_plot_dist_args(subparser): Add common arguments for distribution graphs
    '''
    # dist(): Required argument
    subparser.add_argument("-i", "--df", help="Input dataframe file path", type=str, required=True)
    subparser.add_argument("-x", "--x", type=str, help="X-axis column name", required=True)

    # File output
    subparser.add_argument("-o", "--dir", type=str, help="Output directory", default='./out')
    subparser.add_argument("-f", "--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_dist.png')

    # Optional core arguments
    subparser.add_argument("-c", "--cols", type=str, help="Color column name for grouping")
    subparser.add_argument("-co", "--cols_ord", nargs="+", help="Custom order for color column values")
    subparser.add_argument("-ce", "--cols_exclude", nargs="+", help="Color column values to exclude")

    # Plot customization
    subparser.add_argument("-b", "--bins", type=int, default=40, help="Number of bins for histogram")
    subparser.add_argument("-LL", "--log10_low", type=int, default=0, help="Log10 scale lower bound (e.g., 1 = 10^1 = 10)")
    subparser.add_argument("-pc", "--palette_or_cmap", type=str, default="colorblind", help="Seaborn color palette or matplotlib colormap")
    subparser.add_argument("-ec", "--edgecol", type=str, default="black", help="Edge color of histogram bars")
    subparser.add_argument("-lw", "--lw", type=int, default=1, help="Line width for edges")
    subparser.add_argument("-ht", "--ht", type=float, default=1.5, help="Height of the plot")
    subparser.add_argument("-asp", "--asp", type=int, default=5, help="Aspect ratio of the plot")
    subparser.add_argument("-tp", "--tp", type=float, default=0.8, help="Top padding space")
    subparser.add_argument("-hs", "--hs", type=int, default=0, help="Horizontal spacing between plots (if faceted)")
    subparser.add_argument("-despine", "--despine", action="store_true", help="Remove plot spines (despine)", default=False)

    # Figure appearance
    subparser.add_argument("-fs", "--figsize", type=parse_tuple_int, default=(5,5), help="Figure size formatted as 'width,height'")
    subparser.add_argument("-t", "--title", type=str, default="", help="Plot title text")
    subparser.add_argument("-ts", "--title_size", type=int, default=18, help="Plot title font size")
    subparser.add_argument("-tw", "--title_weight", type=str, default="bold", help="Plot title font weight (e.g., bold, normal)")
    subparser.add_argument("-tf", "--title_font", type=str, default="Arial", help="Font family for the plot title")

    # X-axis
    subparser.add_argument("-xa", "--x_axis", type=str, default="", help="Label for the X-axis")
    subparser.add_argument("-xs", "--x_axis_size", type=int, default=12, help="Font size for X-axis label")
    subparser.add_argument("-xw", "--x_axis_weight", type=str, default="bold", help="Font weight for X-axis label")
    subparser.add_argument("-xf", "--x_axis_font", type=str, default="Arial", help="Font family for X-axis label")
    subparser.add_argument("-xsc", "--x_axis_scale", type=str, default="linear", help="X-axis scale (e.g., linear, log)")
    subparser.add_argument("-xd", "--x_axis_dims", type=parse_tuple_float, default=(0, 0), help="X-axis range as tuple: start,end")
    subparser.add_argument("-xp", "--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("-xts", "--x_ticks_size", type=int, default=9, help="Font size for X-axis tick labels")
    subparser.add_argument("-xtr", "--x_ticks_rot", type=int, default=0, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("-xtf", "--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")
    subparser.add_argument("-xt", "--x_ticks", nargs="+", help="Explicit tick values for X-axis")
    # Y-axis
    subparser.add_argument("-ya", "--y_axis", type=str, default="", help="Label for the Y-axis")
    subparser.add_argument("-ys", "--y_axis_size", type=int, default=12, help="Font size for Y-axis label")
    subparser.add_argument("-yw", "--y_axis_weight", type=str, default="bold", help="Font weight for Y-axis label")
    subparser.add_argument("-yf", "--y_axis_font", type=str, default="Arial", help="Font family for Y-axis label")
    subparser.add_argument("-ysc", "--y_axis_scale", type=str, default="linear", help="Y-axis scale (e.g., linear, log)")
    subparser.add_argument("-yd", "--y_axis_dims", type=parse_tuple_float, default=(0, 0), help="Y-axis range as tuple: start,end")
    subparser.add_argument("-yp", "--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("-yts", "--y_ticks_size", type=int, default=9, help="Font size for Y-axis tick labels")
    subparser.add_argument("-ytr", "--y_ticks_rot", type=int, default=0, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("-ytf", "--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")
    subparser.add_argument("-yt", "--y_ticks", nargs="+", help="Explicit tick values for Y-axis")

    # Legend
    subparser.add_argument("-lt", "--legend_title", type=str, default="", help="Title text for the legend")
    subparser.add_argument("-lts", "--legend_title_size", type=int, default=12, help="Font size of the legend title")
    subparser.add_argument("-ls", "--legend_size", type=int, default=9, help="Font size for legend items")
    subparser.add_argument("-lba", "--legend_bbox_to_anchor", type=parse_tuple_float, default=(1, 1), help="Legend bbox anchor position")
    subparser.add_argument("-ll", "--legend_loc", type=str, default="upper left", help="Legend location on the plot")
    subparser.add_argument("-li", "--legend_items", type=parse_tuple_int, default=(0, 0), help="Tuple for legend layout items")
    subparser.add_argument("-ln", "--legend_ncol", type=int, default=1, help="Number of columns in the legend")

    # Final display
    subparser.add_argument("-d", "--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
    subparser.add_argument("-s", "--show", action="store_true", help="Show the plot in an interactive window", default=False)
    subparser.add_argument("-sc", "--space_capitalize", action="store_true", help="Capitalize and space legend/label values", default=False)

def add_common_plot_heat_args(subparser, fastq_parser=False, stat_parser=False):
    '''
    add_common_plot_heat_args(subparser): Add common arguments for heatmap graphs
    '''
    # Required arguments
    if stat_parser == False:
        subparser.add_argument("-i", "--df", help="Input dataframe file path", type=str, required=True)

    # Optional arguments
    if fastq_parser == False and stat_parser == False:
        subparser.add_argument("-x", "--x", type=str, help="X-axis column name to pivot tidy-formatted dataframe into matrix format")
        subparser.add_argument("-y", "--y", type=str, help="Y-axis column name to pivot tidy-formatted dataframe into matrix format")
        subparser.add_argument("-vr", "--vars", type=str, help="Variable column name to split tidy-formatted dataframe into a dictionary of pivoted dataframes")
        subparser.add_argument("-vs", "--vals", type=str, help="Value column name to populate pivoted dataframes")
    if fastq_parser == False:    
        subparser.add_argument("-vd", "--vals_dims", type=parse_tuple_float, help="Value column limits formatted as 'vmin,vmax'")
    if fastq_parser == True:
        subparser.add_argument("-cc", "--cond_col", type=str, help="Condition column name", required=True)
        subparser.add_argument("-c", "--cond", type=str, help="Condition value for filtering", required=True)
        subparser.add_argument("-FC", "--FC", type=str, help="Fold change column name (values within heatmap after log2 transformation)", required=True)
        subparser.add_argument("-wp", "--wt_prot", type=str, help="WT protein sequence", required=True)
        subparser.add_argument("-wr", "--wt_res", type=int, help="WT protein sequence residue start number", required=True)
        subparser.add_argument("-co", "--cutoff", type=float, help="Comparison count mean cutoff for masking low-abundance values", default=argparse.SUPPRESS)
        subparser.add_argument("-aa", "--aa", type=str, help="AA saturation mutagenesis. The 'aa' option [default] makes all amino acid substitutions ('aa_subs'), +1 amino acid insertions ('aa_ins'), and -1 amino acid deletions ('aa_dels').", choices=['aa', 'aa_subs', 'aa_ins', 'aa_dels'], default='aa')
        subparser.add_argument("-zc", "--z_col", type=str, help="Column name for Z-score normalization (Default: None)", default=argparse.SUPPRESS)
        subparser.add_argument("-zv", "--z_var", type=str, help="Column value for Z-score normalization (Default: None)", default=argparse.SUPPRESS)
        subparser.add_argument("-l", "--label", type=str, help="Label column name (Default: 'Edit'). Can't be None.", default='Edit')
    
    if stat_parser == False:
        subparser.add_argument("-o", "--dir", type=str, help="Output directory path", default='./out')
        subparser.add_argument("-f", "--file", type=str, help="Output filename", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_heat.png')
    subparser.add_argument("-ec", "--edgecol", type=str, default="black", help="Color of cell edges")
    subparser.add_argument("-lw", "--lw", type=int, default=1, help="Line width for cell borders")

    if fastq_parser == False:
        subparser.add_argument("-a", "--annot", action="store_true", help="Display cell values as annotations", default=False)
    subparser.add_argument("-center", "--center", type=float, default=0, help="Center value for colormap (Default: 0)")
    subparser.add_argument("-cm", "--cmap", type=str, default="Reds", help="Matplotlib colormap to use for heatmap")
    if fastq_parser == True:
        subparser.add_argument("-cmWT", "--cmap_WT", type=str, default="forestgreen", help="Color for WT values in colorbar (Default: forestgreen)")
        subparser.add_argument("-cmnWT", "--cmap_not_WT", type=str, default="lightgray", help="Color for non-WT values in colorbar (Default: lightgray)")
    subparser.add_argument("-sq", "--sq", action="store_true", help="Use square aspect ratio for cells", default=False)
    subparser.add_argument("-cb", "--cbar", action="store_true", help="Display colorbar", default=False)
    if stat_parser == False:
        subparser.add_argument("-cl", "--cbar_label", type=str, default=argparse.SUPPRESS, help="Colorbar label")
    subparser.add_argument("-cls", "--cbar_label_size", type=int, default=argparse.SUPPRESS, help="Font size for colorbar label")
    subparser.add_argument("-clw", "--cbar_label_weight", type=str, default="bold", help="Font weight for colorbar label (Default: bold)", choices=['bold', 'normal', 'heavy'])
    subparser.add_argument("-cts", "--cbar_tick_size", type=int, default=argparse.SUPPRESS, help="Font size for colorbar ticks")
    subparser.add_argument("-cs", "--cbar_shrink", type=float, default=argparse.SUPPRESS, help="Shrink factor for colorbar")
    subparser.add_argument("-ca", "--cbar_aspect", type=int, default=argparse.SUPPRESS, help="Aspect ratio for colorbar")
    subparser.add_argument("-cp", "--cbar_pad", type=float, default=argparse.SUPPRESS, help="Padding for colorbar")
    subparser.add_argument("-cor", "--cbar_orientation", type=str, default=argparse.SUPPRESS, help="Orientation of colorbar (Default: 'vertical')", choices=['vertical', 'horizontal'])

    # Title and size
    subparser.add_argument("-t", "--title", type=str, default="", help="Plot title")
    subparser.add_argument("-ts", "--title_size", type=int, default=18, help="Font size of the title")
    subparser.add_argument("-tw", "--title_weight", type=str, default="bold", help="Font weight of the title (e.g., bold, normal)")
    subparser.add_argument("-tf", "--title_font", type=str, default="Arial", help="Font family for the title")
    subparser.add_argument("-fs", "--figsize", type=parse_tuple_int, default=(5,5), help="Figure size formatted as 'width,height'")

    # X-axis
    subparser.add_argument("-xa", "--x_axis", type=str, default="", help="X-axis label")
    subparser.add_argument("-xas", "--x_axis_size", type=int, default=12, help="Font size for X-axis label")
    subparser.add_argument("-xaw", "--x_axis_weight", type=str, default="bold", help="Font weight for X-axis label")
    subparser.add_argument("-xaf", "--x_axis_font", type=str, default="Arial", help="Font family for X-axis label")
    subparser.add_argument("-xap", "--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("-xts", "--x_ticks_size", type=int, default=9, help="Font size for X-axis tick labels")
    subparser.add_argument("-xtr", "--x_ticks_rot", type=int, default=0, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("-xtf", "--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")

    # Y-axis
    subparser.add_argument("-ya", "--y_axis", type=str, default="", help="Y-axis label")
    subparser.add_argument("-yas", "--y_axis_size", type=int, default=12, help="Font size for Y-axis label")
    subparser.add_argument("-yaw", "--y_axis_weight", type=str, default="bold", help="Font weight for Y-axis label")
    subparser.add_argument("-yaf", "--y_axis_font", type=str, default="Arial", help="Font family for Y-axis label")
    subparser.add_argument("-yap", "--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("-yts", "--y_ticks_size", type=int, default=9, help="Font size for Y-axis tick labels")
    subparser.add_argument("-ytr", "--y_ticks_rot", type=int, default=0, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("-ytf", "--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")

    # Final display
    subparser.add_argument("-d", "--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
    subparser.add_argument("-s","--show", action="store_true", help="Show the plot in an interactive window", default=False)
    if fastq_parser == False:
        subparser.add_argument("--space_capitalize", action="store_true", help="Capitalize and space labels/legend values", default=False)

def add_common_plot_stack_args(subparser, fastq_parser=False):
    '''
    add_common_plot_stack_args(subparser): Add common arguments for stacked bar plot
    '''
    # Required arguments
    subparser.add_argument("-i", "--df", type=str, help="Input dataframe file path", required=True)
    subparser.add_argument("-x", "--x", type=str, help="X-axis column name")
    subparser.add_argument("-y", "--y", type=str, help="Y-axis column name")
    subparser.add_argument("-c", "--cols", type=str, help="Color column name for stacking")

    # Optional parameters
    subparser.add_argument("-o", "--dir", type=str, help="Output directory path", default='./out')
    subparser.add_argument("-f", "--file", type=str, help="Output filename", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_stack.png')
    
    subparser.add_argument("-cg", "--cutoff_group", type=str, default=argparse.SUPPRESS, help="Column name to group by when applying cutoff")
    subparser.add_argument("-cv", "--cutoff_value", type=float, default=0, help="Y-axis values needs be greater than (e.g. 0)")
    subparser.add_argument("-cr", "--cutoff_remove", dest="cutoff_keep",action="store_false", help="Remove values below cutoff", default=True)
    subparser.add_argument("-co", "--cols_ord", nargs="+", help="Order of values in the color column")
    subparser.add_argument("-xo", "--x_ord", nargs="+", help="Custom order of X-axis categories")
    subparser.add_argument("-pc", "--palette_or_cmap", type=str, default="Set2", help="Seaborn palette or Matplotlib colormap for stacked bars")
    subparser.add_argument("-ec", "--errcap", type=int, default=4, help="Width of error bar caps")
    subparser.add_argument("-v", "--vertical", action="store_true", help="Stack bars vertically (default True)", default=False)
    if fastq_parser==True:
        subparser.add_argument("-PDB_pt", "--PDB_pt", type=str, default=argparse.SUPPRESS, help="PDB ID (if saved to ~/.config/edms/PDB) or file path for PDB structure file. See edms.dat.pdb.retrieve() or edms uniprot retrieve -h for more information")

    # Figure & layout
    subparser.add_argument("-fs", "--figsize", type=parse_tuple_int, default=(5,5), help="Figure size formatted as 'width,height'")
    subparser.add_argument("-t", "--title", type=str, default="", help="Plot title")
    subparser.add_argument("-ts", "--title_size", type=int, default=18, help="Font size of the title")
    subparser.add_argument("-tw", "--title_weight", type=str, default="bold", help="Font weight of the title (e.g., bold, normal)")
    subparser.add_argument("-tf", "--title_font", type=str, default="Arial", help="Font family for the title")

    # X-axis formatting
    subparser.add_argument("-xa", "--x_axis", type=str, default="", help="X-axis label")
    subparser.add_argument("-xas", "--x_axis_size", type=int, default=12, help="Font size for X-axis label")
    subparser.add_argument("-xaw", "--x_axis_weight", type=str, default="bold", help="Font weight for X-axis label")
    subparser.add_argument("-xaf", "--x_axis_font", type=str, default="Arial", help="Font family for X-axis label")
    subparser.add_argument("-xap", "--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("-xts", "--x_ticks_size", type=int, default=9, help="Font size for X-axis tick labels")
    subparser.add_argument("-xtr", "--x_ticks_rot", type=int, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("-xtf", "--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")

    # Y-axis formatting
    subparser.add_argument("-ya", "--y_axis", type=str, default="", help="Y-axis label")
    subparser.add_argument("-yas", "--y_axis_size", type=int, default=12, help="Font size for Y-axis label")
    subparser.add_argument("-yaw", "--y_axis_weight", type=str, default="bold", help="Font weight for Y-axis label")
    subparser.add_argument("-yaf", "--y_axis_font", type=str, default="Arial", help="Font family for Y-axis label")
    subparser.add_argument("-yad", "--y_axis_dims", type=parse_tuple_float, default=(0,0), help="Y-axis range as tuple: start,end")
    subparser.add_argument("-yap", "--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("-yts", "--y_ticks_size", type=int, default=9, help="Font size for Y-axis tick labels")
    subparser.add_argument("-ytr", "--y_ticks_rot", type=int, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("-ytf", "--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")

    # Legend options
    subparser.add_argument("-lt", "--legend_title", type=str, default="", help="Legend title text")
    subparser.add_argument("-lts", "--legend_title_size", type=int, default=12, help="Font size of the legend title")
    subparser.add_argument("-ls", "--legend_size", type=int, default=12, help="Font size for legend items")
    subparser.add_argument("-lba", "--legend_bbox_to_anchor", type=parse_tuple_float, default=(1, 1), help="Anchor position for the legend bounding box")
    subparser.add_argument("-ll", "--legend_loc", type=str, default="upper left", help="Legend location on the plot")
    subparser.add_argument("-ln", "--legend_ncol", type=int, default=1, help="Number of columns in the legend")
    subparser.add_argument("-lcs", "--legend_columnspacing", type=int, default=argparse.SUPPRESS, help="Space between columns in legend; only for html plots")
    subparser.add_argument("-lhtp", "--legend_handletextpad", type=float, default=argparse.SUPPRESS, help="Space between marker and text in legend; only for html plots")
    subparser.add_argument("-lls", "--legend_labelspacing", type=float, default=argparse.SUPPRESS, help="Vertical space between entries in legend; only for html plots")
    subparser.add_argument("-lbp", "--legend_borderpad", type=float, default=argparse.SUPPRESS, help="Padding inside legend box; only for html plots")
    subparser.add_argument("-lhl", "--legend_handlelength", type=float, default=argparse.SUPPRESS, help="Marker length in legend; only for html plots")
    subparser.add_argument("-lshm", "--legend_size_html_multiplier", type=float, default=argparse.SUPPRESS, help="Legend size multiplier for html plots")

    # Display and formatting
    subparser.add_argument("-d","--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
    subparser.add_argument("-s","--show", action="store_true", help="Show the plot in an interactive window", default=False)
    subparser.add_argument("-sc","--space_capitalize", action="store_true", help="Capitalize and space legend/label values", default=False)

def add_common_plot_vol_args(subparser, fastq_parser=False):
    '''
    add_common_plot_vol_args(subparser): Add common arguments for volcano plot
    '''
    # Required arguments
    subparser.add_argument("-i", "--df", type=str, help="Input dataframe file path from edms stat compare", required=True)
    subparser.add_argument("-FC", "--FC", type=str, help="Fold change column name (X-axis)", required=True)
    subparser.add_argument("-pval", "--pval", type=str, help="P-value column name (Y-axis)", required=True)

    # Optional data columns
    subparser.add_argument("-stys", "--stys", type=str, help="Style column name for custom markers")
    subparser.add_argument("-size", "--size", type=str, help="Column name used to scale point sizes (Default: -log10('pval'); specify 'false' for no size)")
    subparser.add_argument("-sd", "--size_dims", type=parse_tuple_float, help="Size range for points formatted as min,max")
    subparser.add_argument("-zc", "--z_col", type=str, help="Column name for Z-score normalization (Default: None)", default=argparse.SUPPRESS)
    subparser.add_argument("-zv", "--z_var", type=str, help="Column value for Z-score normalization (Default: None)", default=argparse.SUPPRESS)    
    subparser.add_argument("-l", "--label", type=str, help="Column containing text labels for points")
    if fastq_parser==False:
        subparser.add_argument("-so","--stys_order", type=str, nargs="+", help="Style column values order")
        subparser.add_argument("-mo", "--mark_order", type=str, nargs="+", help="Markers order for style column values order")
    
    # Additional annotation data sources
    if fastq_parser==True:
        add_common_fastq_label_args(subparser)
        
    # Thresholds
    subparser.add_argument("-ft","--FC_threshold", type=float, default=1, help="Fold change threshold for significance")
    subparser.add_argument("-pt","--pval_threshold", type=float, default=1, help="P-value threshold for significance")

    # Output
    subparser.add_argument("-o", "--dir", type=str, help="Output directory path", default='./out')
    subparser.add_argument("-f", "--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_plot_vol.png')

    # Aesthetics
    subparser.add_argument("-c","--color", type=str, default="lightgray", help="Color for non-significant points")
    subparser.add_argument("-a","--alpha", type=float, default=0.5, help="Transparency for non-significant points")
    subparser.add_argument("-ec","--edgecol", type=str, default="black", help="Edge color of points")
    subparser.add_argument("-v","--vertical", action="store_true", help="Use vertical layout for plot", default=False)

    # Figure setup
    subparser.add_argument("-fs","--figsize", type=parse_tuple_int, default=(5,5), help="Figure size formatted as 'width,height'")
    subparser.add_argument("-t","--title", type=str, default="", help="Plot title")
    subparser.add_argument("-ts","--title_size", type=int, default=18, help="Font size for plot title")
    subparser.add_argument("-tw","--title_weight", type=str, default="bold", help="Font weight for plot title (e.g., bold, normal)")
    subparser.add_argument("-tf","--title_font", type=str, default="Arial", help="Font family for plot title")
    # X-axis settings
    subparser.add_argument("-xa","--x_axis", type=str, default="", help="Label for the X-axis")
    subparser.add_argument("-xas","--x_axis_size", type=int, default=12, help="Font size for X-axis label")
    subparser.add_argument("-xaw","--x_axis_weight", type=str, default="bold", help="Font weight for X-axis label")
    subparser.add_argument("-xaf","--x_axis_font", type=str, default="Arial", help="Font family for X-axis label")
    subparser.add_argument("-xad","--x_axis_dims", type=parse_tuple_float, default=(0, 0), help="X-axis range formatted as min,max")
    subparser.add_argument("-xap","--x_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for X-axis label")
    subparser.add_argument("-xts","--x_ticks_size", type=int, default=9, help="Font size for X-axis tick labels")
    subparser.add_argument("-xtr","--x_ticks_rot", type=int, default=0, help="Rotation angle for X-axis tick labels")
    subparser.add_argument("-xtf","--x_ticks_font", type=str, default="Arial", help="Font family for X-axis tick labels")
    subparser.add_argument("-xt","--x_ticks", nargs="+", help="Custom tick values for X-axis")

    # Y-axis settings
    subparser.add_argument("-ya","--y_axis", type=str, default="", help="Label for the Y-axis")
    subparser.add_argument("-yas","--y_axis_size", type=int, default=12, help="Font size for Y-axis label")
    subparser.add_argument("-yaw","--y_axis_weight", type=str, default="bold", help="Font weight for Y-axis label")
    subparser.add_argument("-yaf","--y_axis_font", type=str, default="Arial", help="Font family for Y-axis label")
    subparser.add_argument("-yad","--y_axis_dims", type=parse_tuple_float, default=(0, 0), help="Y-axis range formatted as min,max")
    subparser.add_argument("-yap","--y_axis_pad", type=int, default=argparse.SUPPRESS, help="Padding for Y-axis label")
    subparser.add_argument("-yts","--y_ticks_size", type=int, default=9, help="Font size for Y-axis tick labels")
    subparser.add_argument("-ytr","--y_ticks_rot", type=int, default=0, help="Rotation angle for Y-axis tick labels")
    subparser.add_argument("-ytf","--y_ticks_font", type=str, default="Arial", help="Font family for Y-axis tick labels")
    subparser.add_argument("-yt","--y_ticks", nargs="+", help="Custom tick values for Y-axis")

    # Legend
    subparser.add_argument("-lt","--legend_title", type=str, default="", help="Title for the legend")
    subparser.add_argument("-lts","--legend_title_size", type=int, default=12, help="Font size for legend title")
    subparser.add_argument("-ls","--legend_size", type=int, default=9, help="Font size for legend items")
    subparser.add_argument("-lba","--legend_bbox_to_anchor", type=parse_tuple_float, default=(1, 1), help="Bounding box anchor for legend")
    subparser.add_argument("-ll","--legend_loc", type=str, default="upper left", help="Legend location on the plot")
    subparser.add_argument("-ln","--legend_ncol", type=int, default=1, help="Number of columns in the legend")
    subparser.add_argument("-lcs",'--legend_columnspacing', type=int, default=argparse.SUPPRESS, help='space between columns in legend; only for html plots')
    subparser.add_argument("-lht",'--legend_handletextpad', type=float, default=argparse.SUPPRESS, help='space between marker and text in legend; only for html plots')
    subparser.add_argument("-lls",'--legend_labelspacing', type=float, default=argparse.SUPPRESS, help='vertical space between entries in legend; only for html plots')
    subparser.add_argument("-lbp",'--legend_borderpad', type=float, default=argparse.SUPPRESS, help='padding inside legend box; only for html plots')
    subparser.add_argument("-lhl",'--legend_handlelength', type=float, default=argparse.SUPPRESS, help='marker length in legend; only for html plots')
    subparser.add_argument("-lshm",'--legend_size_html_multiplier', type=float, default=argparse.SUPPRESS, help='legend size multiplier for html plots')
    
    # Boolean switches
    subparser.add_argument("-ddl","--dont_display_legend", action="store_false", help="Don't display legend on plot", default=True)
    subparser.add_argument("-dl","--display_labels", type=str, nargs="+", help="Display labels for values if label column specified (Options: 'FC & p-value', 'FC', 'p-value', 'NS', 'all', or ['label1', 'label2', ..., 'labeln'])", default=["FC & p-value"])
    subparser.add_argument("-dda","--dont_display_axis", dest='display_axis', action="store_false", default=True, help="Display x- and y-axis lines (Default: True)")
    subparser.add_argument("-dli","--display_lines", action="store_true", help="Display lines for threshold (Default: False)", default=False)
    subparser.add_argument("-d","--dpi", type=int, help="Figure dpi (Default: 600 for non-HTML, 150 for HTML)", default=0)
    subparser.add_argument("-s","--show", action="store_true", help="Show the plot in an interactive window", default=False)
    subparser.add_argument("-sc","--space_capitalize", action="store_true", help="Capitalize and space labels/legend items", default=False)

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
    - odds_ratio() [or]: computes odds ratio relative to a specified condition (OR = (A/B)/(C/D))
    - zscore(): Z-score `val` within each `cond_col` using stats computed from rows where `var_col == var`.
    '''
    parser_stat = subparsers.add_parser("stat", help="Statistics", description="Statistics", formatter_class=formatter_class)
    subparsers_stat = parser_stat.add_subparsers()
    
    # describe(): returns descriptive statistics for numerical columns in a DataFrame
    parser_stat_describe = subparsers_stat.add_parser("describe", help="Compute descriptive statistics", description="Compute descriptive statistics", formatter_class=formatter_class)

    parser_stat_describe.add_argument("-i", "--df", type=str, help="Input file path", required=True)

    parser_stat_describe.add_argument("-o", "--dir", type=str, help="Output directory (Default: ../out)",default='../out')
    parser_stat_describe.add_argument("-f", "--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_descriptive.csv')
    
    parser_stat_describe.add_argument("-c","--cols", nargs="+", help="List of numerical columns to describe")
    parser_stat_describe.add_argument("-g","--group", type=str, help="Column name to group by")
    
    parser_stat_describe.set_defaults(func=st.describe)

    # difference(): computes the appropriate statistical test(s) and returns the p-value(s)
    parser_stat_difference = subparsers_stat.add_parser("difference", help="Compute statistical difference between groups", description="Compute statistical difference between groups", formatter_class=formatter_class)

    parser_stat_difference.add_argument("-i", "--df", type=str, help="Input file path",required=True)
    parser_stat_difference.add_argument("-d","--data_col", type=str, help="Name of column containing numerical data",required=True)
    parser_stat_difference.add_argument("-c","--compare_col", type=str, help="Name of column used for grouping/comparisons",required=True)
    parser_stat_difference.add_argument("-p","--compare", nargs="+", help="List of groups to compare (e.g. A B)",required=True)

    parser_stat_difference.add_argument("-o", "--dir", type=str, help="Output directory (Default: ../out)",default='../out')
    parser_stat_difference.add_argument("-f", "--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_difference.csv')

    parser_stat_difference.add_argument("-s","--same", action="store_true", help="Same subjects (paired test)")
    parser_stat_difference.add_argument("-r","--para", action="store_true", help="Use parametric test (Default: True)")
    parser_stat_difference.add_argument("-a","--alpha", type=float, default=0.05, help="Significance level (Default: 0.05)")
    parser_stat_difference.add_argument("-w","--within_cols", nargs="+", help="Columns for repeated measures (used if same=True and para=True)")
    parser_stat_difference.add_argument("-m","--method", type=str, default="holm", help="Correction method for multiple comparisons")
    parser_stat_difference.set_defaults(func=st.difference)

    # correlation(): returns a correlation matrix
    parser_stat_correlation = subparsers_stat.add_parser("correlation", help="Compute correlation matrix", description="Compute correlation matrix", formatter_class=formatter_class)

    parser_stat_correlation.add_argument("-i", "--df", type=str, help="Input file path",required=True)
    parser_stat_correlation.add_argument("-v","--var_cols", nargs="+", help="List of 2 variable columns for tidy format")
    parser_stat_correlation.add_argument("-V","--value_cols", nargs="+", help="List of numerical columns to correlate")
    parser_stat_correlation.add_argument("-m","--method", type=str, default="pearson", choices=["pearson", "spearman", "kendall"],
                                         help="Correlation method to use (Default: pearson)")
    parser_stat_correlation.add_argument("-n","--numeric_only", action="store_true", help="Only use numeric columns (Default: True)")
    parser_stat_correlation.add_argument("-N","--no_plot", dest="plot", action="store_false", help="Don't generate correlation matrix plot", default=True)
    parser_stat_correlation.add_argument("-o", "--dir", type=str, help="Output directory (Default: ../out)",default='../out')
    parser_stat_correlation.add_argument("-F", "--file_data", type=str, help="Output data file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_correlation.csv')
    parser_stat_correlation.add_argument("-P", "--file_plot", type=str, help="Output plot file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_correlation.pdf')
    add_common_plot_heat_args(parser_stat_correlation, stat_parser=True)

    parser_stat_correlation.set_defaults(func=st.correlation)

    # compare(): computes FC, pval, and log transformations relative to a specified condition
    parser_stat_compare = subparsers_stat.add_parser("compare", help="Compare conditions using FC, p-values, and log transforms", description="Compare conditions using FC, p-values, and log transforms", formatter_class=formatter_class)

    parser_stat_compare.add_argument("-i", "--df", type=str, help="Input file path",required=True)
    parser_stat_compare.add_argument("-s","--sample", type=str, help="Sample column name",required=True)
    parser_stat_compare.add_argument("-c","--cond", type=str, help="Condition column name",required=True)
    parser_stat_compare.add_argument("-C","--cond_comp", type=str, help="Condition to compare against",required=True)
    parser_stat_compare.add_argument("-V","--var", type=str, help="Variable column name",required=True)
    parser_stat_compare.add_argument("-n","--count", type=str, help="Count column name",required=True)
    parser_stat_compare.add_argument("-p","--pseudocount", type=int, default=1, help="Pseudocount to avoid log(0) or divide-by-zero errors")
    parser_stat_compare.add_argument("-r","--replicate", type=str, help="Replicate column name (use pairwise comparisons instead of aggregate)", default=argparse.SUPPRESS)
    parser_stat_compare.add_argument("-a","--alternative", type=str, default="two-sided", choices=["two-sided", "less", "greater"], help="Alternative hypothesis for Fisher's exact test (Default: two-sided)")
    parser_stat_compare.add_argument("-o", "--dir", type=str, help="Output directory (Default: ../out)",default='../out')
    parser_stat_compare.add_argument("-f", "--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_compare.csv')
    parser_stat_compare.add_argument("-v","--verbose", action="store_true", help="Print progress to console", default=False)

    parser_stat_compare.set_defaults(func=st.compare)

    # odds_ratio() [or]: computes odds ratio relative to a specified condition (OR = (A/B)/(C/D))
    parser_stat_odds_ratio = subparsers_stat.add_parser("or", help="Computes odds ratios relative to a specified condition & variable (e.g., unedited & WT)", description="Compute odds ratio between conditions", formatter_class=formatter_class)

    parser_stat_odds_ratio.add_argument("-i", "--df", type=str, help="Input file path",required=True)
    parser_stat_odds_ratio.add_argument("-c","--cond", type=str, help="Condition column name",required=True)
    parser_stat_odds_ratio.add_argument("-C","--cond_comp", type=str, help="Condition for comparison group",required=True)
    parser_stat_odds_ratio.add_argument("-V","--var", type=str, help="Variable column name",required=True)
    parser_stat_odds_ratio.add_argument("-W","--var_comp", type=str, help="Variable name for comparison (e.g., WT)",required=True)
    parser_stat_odds_ratio.add_argument("-n","--count", type=str, help="Count column name",required=True)
    
    parser_stat_odds_ratio.add_argument("-p","--pseudocount", type=int, default=1, help="Pseudocount to avoid /0 (Default: 1)")
    parser_stat_odds_ratio.add_argument("-a","--alternative", type=str, default="two-sided", choices=["two-sided", "less", "greater"], help="Alternative hypothesis for Fisher's exact test (Default: two-sided)")
    parser_stat_odds_ratio.add_argument("-o", "--dir", type=str, help="Output directory (Default: ../out)",default='../out')
    parser_stat_odds_ratio.add_argument("-f", "--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_odds_ratio.csv')
    parser_stat_odds_ratio.add_argument("-v","--verbose", action="store_true", help="Print progress to console", default=False)

    parser_stat_odds_ratio.set_defaults(func=st.odds_ratio)

    # zscore(): Z-score `val` within each `cond_col` using stats computed from rows where `var_col == var`.
    parser_stat_zscore = subparsers_stat.add_parser("zscore", help="Compute Z-scores within conditions using stats from a reference variable", description="Compute Z-scores within conditions using stats from a reference variable", formatter_class=formatter_class)

    parser_stat_zscore.add_argument("-i", "--df", type=str, help="Input file path",required=True)
    parser_stat_zscore.add_argument("-c","--cond_col", type=str, help="Condition column name",required=True)
    parser_stat_zscore.add_argument("-v","--var_col", type=str, help="Variable column name",required=True)
    parser_stat_zscore.add_argument("-V","--var", type=str, help="Variable name for reference stats",required=True) 
    parser_stat_zscore.add_argument("-vl","--val_col", type=str, help="Value column name to Z-score",required=True)

    parser_stat_zscore.add_argument("-oc","--out_col", type=str, help="Name for output Z-score column", default=argparse.SUPPRESS)
    parser_stat_zscore.add_argument("-d","--ddof", type=float, help="Delta degrees of freedom for std (0 = population, 1 = sample [default]).", default=1)
    parser_stat_zscore.add_argument("-o", "--dir", type=str, help="Output directory (Default: ../out)",default='../out')
    parser_stat_zscore.add_argument("-f", "--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_zscore.csv')
    
    parser_stat_zscore.set_defaults(func=st.zscore)

    '''
    edms.gen.io:
    - in_subs() [in]: moves all files with a given suffix into subfolders named after the files (excluding the suffix).
    - out_subs() [out]: recursively moves all files from subdirectories into the parent directory and delete the emptied subdirectories.
    - create_sh() [sh]: creates a shell script with SLURM job submission parameters for Harvard FASRC cluster.
    - combine(): Combine text files matching provided suffixes into a single output file, inserting a header with the original filename before each file's content.
    - split_R1_R2(): split paired reads into new R1 and R2 subdirectories at the parent directory
    - excel_csvs(): exports excel file to .csv files in specified directory  
    '''
    parser_io = subparsers.add_parser("io", help="Input/Output", formatter_class=formatter_class)
    subparsers_io = parser_io.add_subparsers()
    
    # Create subparsers for commands
    parser_io_in_subs = subparsers_io.add_parser("in", help="*No FASRC* Moves all files with a given suffix into subfolders named after the files (excluding the suffix)", description="*No FASRC* Moves all files with a given suffix into subfolders named after the files (excluding the suffix)", formatter_class=formatter_class)
    parser_io_out_subs = subparsers_io.add_parser("out", help="*No FASRC* Delete subdirectories and move their files to the parent directory", description="*No FASRC* Delete subdirectories and move their files to the parent directory", formatter_class=formatter_class)
    parser_io_create_sh = subparsers_io.add_parser("sh", help='Generate SLURM shell script for Harvard FASRC cluster.', description='Generate SLURM shell script for Harvard FASRC cluster.', formatter_class=formatter_class)
    parser_io_combine = subparsers_io.add_parser("combine", help='Combine text files matching provided suffixes into a single output file, inserting a header with the original filename before each file\'s content.', description='Combine text files matching provided suffixes into a single output file, inserting a header with the original filename before each file\'s content.', formatter_class=formatter_class)
    parser_io_split_R1_R2 = subparsers_io.add_parser("split_R1_R2", help='*No FASRC* Split paired reads into new R1 and R2 subdirectories at the parent directory.', description='*No FASRC* Split paired reads into new R1 and R2 subdirectories at the parent directory.', formatter_class=formatter_class)
    parser_io_excel_csvs = subparsers_io.add_parser("excel_csvs", help='Exports excel file to .csv files in specified directory.', description='Exports excel file to .csv files in specified directory.', formatter_class=formatter_class)

    # Add common arguments
    for parser_io_common in [parser_io_in_subs,parser_io_out_subs,parser_io_split_R1_R2]:
        parser_io_common.add_argument("-o", "--dir", help="Path to parent directory", type=str, default='.')
    
    # in_subs() [in] arguments
    parser_io_in_subs.add_argument("-s", "--suf", help="File suffix (e.g., '.txt', '.csv') to filter files.", type=str, required=True)
    parser_io_in_subs.add_argument("-g", "--group_by", help="How to group files into subdirectories (Options: 'basename' [Default] or 'prefix')", type=str, choices=['basename', 'prefix'], default='basename')
    parser_io_in_subs.add_argument("-p", "--prefix_sep", help="Delimiter to use when grouping by prefix", type=str, default=argparse.SUPPRESS)
    parser_io_in_subs.add_argument("-l", "--prefix_len", help="Number of characters to use as prefix.", type=int, default=argparse.SUPPRESS)
    
    # create_sh() [sh]: creates a shell script with SLURM job submission parameters for Harvard FASRC cluster.
    parser_io_create_sh.add_argument('-o', '--dir', type=str, help='Directory to save the shell script.', default='.')
    parser_io_create_sh.add_argument('-f', '--file', type=str, help='Name of the shell script file to create.',default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.sh')
    parser_io_create_sh.add_argument('-c', '--cores', type=int, default=1, help='Number of CPU cores to request.')
    parser_io_create_sh.add_argument('-p', '--partition', type=str, default='serial_requeue', help='SLURM partition to use.')
    parser_io_create_sh.add_argument('-t', '--time', type=str, default='0-00:10', help='Job run time in D-HH:MM format.')
    parser_io_create_sh.add_argument('-m', '--mem', type=int, default=1000, help='Memory in MB.')
    parser_io_create_sh.add_argument('-e', '--email', type=str, default=None, help='Notification email address.')
    parser_io_create_sh.add_argument('-y', '--python', type=str, default='python/3.12.5-fasrc01', help='Python module to load.')
    parser_io_create_sh.add_argument('-n', '--env', type=str, default='edms', help='Conda environment to activate.')

    # combine() arguments
    parser_io_combine.add_argument('-i', '--in_dir', type=str, help='Directory to search for input files.', required=True)
    
    parser_io_combine.add_argument('-o', '--out_dir', type=str, help='Directory to write the combined file (Default: working directory).', default='.')
    parser_io_combine.add_argument('-f', '--out_file', type=str, help='Output filename (Default: $date_time$_combined.txt file).', default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_combined.txt')
    parser_io_combine.add_argument('-s', '--suffixes', type=str, nargs='+', help='Iterable of suffixes to match (e.g. .txt .log .out .err).', default=['.txt', '.log', '.out', '.err'])
    parser_io_combine.add_argument('-r', '--recursive', action='store_true', help='If set, search subdirectories recursively.', default=False)
    parser_io_combine.add_argument('-l', '--full_path', action='store_true', help='If set, include full path in header; otherwise only filename.', default=False)
    parser_io_combine.add_argument('-e', '--encoding', type=str, help='Text encoding to use when reading/writing files.', default='utf-8')
    
    # excel_csvs(): exports excel file to .csv files in specified directory 
    parser_io_excel_csvs.add_argument('-p', '--pt', type=str, help='Excel file path', required=True)
    parser_io_excel_csvs.add_argument('-o', '--dir', type=str, help='Output directory path (Default: same directory as excel file name).',default='')

    # Call command functions
    parser_io_in_subs.set_defaults(func=io.in_subs)
    parser_io_out_subs.set_defaults(func=io.out_subs)
    parser_io_create_sh.set_defaults(func=io.create_sh)
    parser_io_combine.set_defaults(func=io.combine)
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
        parser_com_common.add_argument("-p", "--pt", help="Path to parent directory", type=str, default='.')
    
    # Smaller_fastq arguments
    parser_com_smaller_fastq.add_argument("-r", "--reads", help="# of reads per fastq file", type=int, default='100000') 
    parser_com_smaller_fastq.add_argument("-s", "--suf", help="Fastq file suffix", type=int, default=".fastq.gz") 
    
    # create_export_var arguments
    parser_com_create_export_var.add_argument("-n", "--name", help="Name of the environment variable (e.g., MYPROJ)", required=True)
    parser_com_create_export_var.add_argument("-p", "--pt", help="Path the variable should point to (e.g., ~/projects/myproj)", required=True)
    parser_com_create_export_var.add_argument("-s", "--shell", choices=["bash", "zsh"], default=argparse.SUPPRESS, help="Shell type)")
    
    # view_export_var arguments
    parser_com_view_export_vars.add_argument("-s", "--shell", choices=["bash", "zsh"], default=argparse.SUPPRESS, help="Shell type")

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
    
    parser_html.add_argument("-o", "--dir", type=str, help="Directory containing HTML files to index", default=".")
    parser_html.add_argument("-f", "--file", type=str, help="Output HTML index file name", default="index.html")
    parser_html.add_argument("-r", "--recursive", action="store_true", help="Recursively search subdirectories for HTML files", default=False)
    parser_html.add_argument("-e", "--exclude", type=str, nargs="+", help="List of filenames to exclude (case insensitive)", default=[])
    parser_html.add_argument("-s", "--sort", type=str, choices=["title", "name", "mtime"], help="Sort HTML files by 'title', 'name', or 'mtime' (modification time)", default="title")
    parser_html.add_argument("-l", "--label", type=str, choices=["title", "stem", "name"], help="Card label source: 'title' (HTML <title>), 'stem' (filename without suffix), or 'name' (full filename)", default="title")
    parser_html.add_argument("-n", "--no_preview", dest="preview", action="store_false", help="Don't include an iframe preview panel in the index", default=True)
    parser_html.add_argument("-g", "--grid_cols", type=int, help="Number of columns in the responsive grid layout", default=3)
    parser_html.add_argument("-i", "--image_types", type=str, nargs="+", help="List of image file extensions to include (e.g. .png .jpg .gif). If not specified, only .html files are included.", default=None)
    parser_html.add_argument("-H", "--preview_height_px", type=int, help="Height of the preview iframe in pixels", default=900)
    parser_html.add_argument("-I", "--icon", type=str, help="Name of the SVG icon file (without .svg) to use as favicon", default="python")
    
    parser_html.set_defaults(func=ht.make_html_index)