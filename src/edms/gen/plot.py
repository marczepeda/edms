''' 
Module: plot.py
Author: Marc Zepeda
Created: 2024-08-05
Description: Plot generation

Usage:
[HTML Tooltips]
- SafeHTMLTooltip: robust HTML tooltip that safely handles missing labels/targets.
- ClickTooltip: persistent HTML tooltip that appears on point click.

[Supporting methods]
- export_mpld3_molstar_html(): Create a standalone HTML file that combines:
    - an mpld3 interactive Matplotlib plot
    - a Mol* viewer panel for a PDB structure
- save_fig(): save static image and optionally interactive HTML or JSON via mpld3.
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
- cat(): creates categorical graphs
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
import json
import mpld3
from pathlib import Path
import shutil
import importlib.resources as pkg_resources

from ..utils import mkdir
from . import io
import edms.resources.molstar as molstar_pkg
import edms.resources.icon as icon_pkg

# HTML Tooltips
class SafeHTMLTooltip(mpld3.plugins.PluginBase):
    """
    MPLD3 plugin: robust HTML tooltip that safely handles missing labels/targets.
    Replaces the use of mpld3.plugins.PointHTMLTooltip to avoid crashes when labels/targets are null.
    """

    JAVASCRIPT = r"""
mpld3.register_plugin("safe_htmltooltip", SafeHtmlTooltipPlugin);
SafeHtmlTooltipPlugin.prototype = Object.create(mpld3.Plugin.prototype);
SafeHtmlTooltipPlugin.prototype.constructor = SafeHtmlTooltipPlugin;
SafeHtmlTooltipPlugin.prototype.requiredProps = ["id"];
SafeHtmlTooltipPlugin.prototype.defaultProps = {
    labels: [],
    targets: [],
    hoffset: 0,
    voffset: 10
};

function SafeHtmlTooltipPlugin(fig, props) {
    mpld3.Plugin.call(this, fig, props);
}

SafeHtmlTooltipPlugin.prototype.draw = function() {
    var obj = mpld3.get_element(this.props.id, this.fig);

    if (!obj) {
        console.warn("SafeHtmlTooltipPlugin: no mpld3 element with id", this.props.id);
        return;
    }

    // Always use arrays, never null
    var labels  = Array.isArray(this.props.labels)  ? this.props.labels  : [];
    var targets = Array.isArray(this.props.targets) ? this.props.targets : [];

    var tooltip = d3.select("body").append("div")
        .attr("class", "mpld3-tooltip")
        .style("position", "absolute")
        .style("z-index", "10")
        .style("visibility", "hidden");

    obj.elements()
        .on("mouseover", function(d, i){
            if (!labels.length || !labels[i]) return;
            tooltip.html(labels[i])
                   .style("visibility", "visible");
        })
        .on("mousemove", function(d, i){
            tooltip
                .style("top",  d3.event.pageY + this.props.voffset + "px")
                .style("left", d3.event.pageX + this.props.hoffset + "px");
        }.bind(this))
        .on("mousedown.callout", function(d, i){
            // Only try to open if we actually have a target URL
            if (!targets.length || !targets[i]) return;
            window.open(targets[i], "_blank");
        })
        .on("mouseout", function(d, i){
            tooltip.style("visibility", "hidden");
        });
};
"""

    def __init__(self, points, labels, targets=None, hoffset=0, voffset=10):
        """
        Parameters
        ----------
        points : matplotlib PathCollection
            The scatter/collection artist to which the tooltips are attached.
        labels : sequence
            HTML strings, one per point.
        targets : sequence, optional
            Optional list of URLs to open on mousedown; may be empty or None.
        hoffset, voffset : int, optional
            Pixel offsets for tooltip position.
        """
        if targets is None:
            targets = []

        # Normalize labels to strings and handle NaNs gracefully
        labels_list = []
        for v in labels:
            if v is None or (isinstance(v, float) and math.isnan(v)):
                labels_list.append("")
            else:
                labels_list.append(str(v))

        self.points = points
        self.labels = labels_list
        self.targets = list(targets)

        self.dict_ = dict(
            type="safe_htmltooltip",
            id=mpld3.plugins.get_id(points),
            labels=self.labels,
            targets=self.targets,
            hoffset=hoffset,
            voffset=voffset,
        )


# Persistent click-to-toggle tooltip plugin for MPLD3
class ClickTooltip(mpld3.plugins.PluginBase):
    """MPLD3 plugin: persistent HTML tooltip that appears on point click.

    Clicking a point shows its label in a tooltip that remains visible
    until another point is clicked. Clicking a point with an empty label
    hides the tooltip.
    """

    JAVASCRIPT = r"""
mpld3.register_plugin("clicktooltip", ClickTooltipPlugin);
ClickTooltipPlugin.prototype = Object.create(mpld3.Plugin.prototype);
ClickTooltipPlugin.prototype.constructor = ClickTooltipPlugin;
ClickTooltipPlugin.prototype.requiredProps = ["id"];
ClickTooltipPlugin.prototype.defaultProps = {
    labels: [],
    hoffset: 0,
    voffset: 10
};

function ClickTooltipPlugin(fig, props) {
    mpld3.Plugin.call(this, fig, props);
}

ClickTooltipPlugin.prototype.draw = function() {
    var obj = mpld3.get_element(this.props.id, this.fig);

    if (!obj) {
        console.warn("ClickTooltip: no mpld3 element with id", this.props.id);
        return;
    }

    var labels = Array.isArray(this.props.labels) ? this.props.labels : [];
    var hoffset = this.props.hoffset || 0;
    var voffset = this.props.voffset || 10;

    // Single tooltip div that we reuse for every click
    var tooltip = d3.select("body").append("div")
        .attr("class", "mpld3-tooltip")
        .style("position", "absolute")
        .style("z-index", "10")
        .style("visibility", "hidden");

    // Click on a point: show/hide tooltip & stop event so the background
    // handler does not immediately hide it.
    obj.elements().on("click", function(d, i) {
        if (d3.event && typeof d3.event.stopPropagation === "function") {
            d3.event.stopPropagation();
        }

        // d3.v3/v5 pattern: `d3.event` holds the mouse event
        var label = labels.length > i ? labels[i] : "";
        if (label == null) label = "";
        label = String(label);

        var hasLabel = label.trim().length > 0;

        if (hasLabel) {
            tooltip.html(label)
                   .style("visibility", "visible")
                   .style("top",  (d3.event.pageY + voffset) + "px")
                   .style("left", (d3.event.pageX + hoffset) + "px");
        } else {
            // Clicking a point with no label hides any previous tooltip
            tooltip.style("visibility", "hidden");
        }
    });

    // Click anywhere else in the figure: hide tooltip
    var canvas = this.fig.canvas;
    if (canvas && typeof canvas.on === "function") {
        canvas.on("click.clicktooltip-bg", function() {
            var evt = d3.event || window.event;
            var target = evt && evt.target ? evt.target : null;

            // Determine if the click target is one of our points
            var nodes;
            if (typeof obj.elements().nodes === "function") {
                // d3 v4/v5 selection
                nodes = obj.elements().nodes();
            } else {
                // d3 v3 selection: first group of elements
                nodes = obj.elements()[0] || [];
            }

            var isPoint = false;
            for (var j = 0; nodes && j < nodes.length; j++) {
                if (nodes[j] === target) {
                    isPoint = true;
                    break;
                }
            }

            if (!isPoint) {
                tooltip.style("visibility", "hidden");
            }
        });
    }
};
"""

    def __init__(self, points, labels, hoffset: int = 0, voffset: int = 10):
        # Normalize labels similarly to SafeHTMLTooltip
        labels_list = []
        for v in labels:
            if v is None or (isinstance(v, float) and math.isnan(v)):
                labels_list.append("")
            else:
                labels_list.append(str(v))

        self.points = points
        self.labels = labels_list
        self.dict_ = dict(
            type="clicktooltip",
            id=mpld3.plugins.get_id(points),
            labels=self.labels,
            hoffset=hoffset,
            voffset=voffset,
        )

# Supporting Methods
def export_mpld3_molstar_html(
    fig: plt.Figure,
    dir: str,
    file: str,
    pdb_id: str | None = None,
    pdb_url: str | None = None,
    icon: str = 'python',
) -> Path:
    """
    export_mpld3_molstar_html(): Create a standalone HTML file that combines:
      - an mpld3 interactive Matplotlib plot
      - a Mol* viewer panel for a PDB structure

    Parameters
    fig (matplotlib.figure.Figure): The Matplotlib figure to convert via mpld3.
    output_html (str or Path): Target HTML file path.
    pdb_id (str, optional): PDB ID to load from PDBe if pdb_url is not given.
    pdb_url (str, optional): Direct URL to a structure file. If given, overrides pdb_id.
    icon (str, optional): html file icon (Default: python logo)

    Returns
    -------
    Path
        Path to the generated HTML file.
    """
    output_html = Path(os.path.join(dir, file))

    # 1) Get mpld3 HTML for the figure (this is just a <div> + <script>, no full document)
    fig_html = mpld3.fig_to_html(fig)

    # 2) Mol* assets (use CDN for standalone HTML)
    molstar_head = f"""
<link rel="stylesheet" type="text/css" href="{file[:-5]}/molstar.css" />
<script src="{file[:-5]}/molstar.js"></script>
"""

    # 3) Set up the structure URL
    # If no explicit URL is given, assume a local PDB file named like '7VOX.pdb'
    if pdb_url is None:
        if pdb_id is None:
            raise ValueError("Either pdb_id or pdb_url must be provided")
        pdb_id_str = str(pdb_id).upper()
        pdb_url = f"{pdb_id_str}.pdb"

    # Infer Mol* format from the file extension (default to 'pdb')
    pdb_ext = Path(pdb_url).suffix.lower().lstrip(".")
    if pdb_ext in {"pdb", "cif", "bcif", "mmcif"}:
        pdb_format = pdb_ext
    else:
        pdb_format = "pdb"

    # Decide whether to load from URL (remote) or embed data (local file)
    url_str = str(pdb_url)
    is_remote = url_str.startswith("http://") or url_str.startswith("https://")

    if is_remote:
        # Keep using URL-based loading; suitable for HTTP(S) servers with CORS enabled
        load_js = f"            viewer.loadStructureFromUrl('{url_str}', '{pdb_format}');"
    else:
        # Embed local file contents directly into the HTML so it works with file:// URLs
        pdb_path = Path(pdb_url)
        if not pdb_path.is_file():
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")
        pdb_text = pdb_path.read_text()
        pdb_js = json.dumps(pdb_text)
        load_js = (
            f"            const pdbData = {pdb_js};\n"
            f"            viewer.loadStructureFromData(pdbData, '{pdb_format}');"
        )

    # 4) Mol* viewer container + initialization script
    molstar_body = rf"""
<div style="display: flex; align-items: flex-start; width: 100%; height: 100%;">
  <div id="mpld3-container" style="width: 50%; min-width: 50%; max-width: 50%; padding-right: 20px; overflow: hidden;">
    {fig_html}
  </div>
  <div id="molstar-container" style="width: 50%; min-width: 50%; max-width: 50%; margin-left: 20px; display: flex; flex-direction: column;">
    <div id="molstar-viewer" style="position: relative; width: 95%; height: 600px;"></div>
  </div>
</div>

<script>
// Keep a global handle for the Mol* viewer so mpld3 callbacks can access it
window.molstarViewer = null;

// Initialize Mol* viewer once the DOM is ready
document.addEventListener('DOMContentLoaded', function () {{
    // 'molstar' global comes from molstar.js
    molstar.Viewer
        .create('molstar-viewer', {{
            layoutIsExpanded: false,
            layoutShowControls: true,
            layoutShowRemoteState: false,
            layoutShowSequence: true,
            layoutShowLog: false,
            layoutShowLeftPanel: false
        }})
        .then(function(viewer) {{
            // Store globally so other scripts (e.g. mpld3 plugins) can talk to Mol*
            window.molstarViewer = viewer;

            // Load structure from embedded data or URL
{load_js}
        }});
}});
</script>
"""

    # 5) Assemble full HTML document
    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8" />
<title>{file[:-5]}</title>
<link rel="icon" type="image/svg+xml" href="{file[:-5]}/{icon}.svg">
{molstar_head}
</head>
<body>
{molstar_body}
</body>
</html>
"""#here

    # 6) Write final file
    output_html.write_text(html, encoding="utf-8")

    return output_html

def save_fig(file: str | None, dir: str | None, fig=None, dpi: int = 0, PDB_pt: str = None, icon: str = 'python') -> None:
    """
    save_fig(): save static image and optionally interactive HTML or JSON via mpld3.

    Parameters:
    file (str | None): output filename (static image by default; `.html` and `.json` trigger interactive mpld3 exports)
    dir (str | None): output directory
    fig: matplotlib Figure object (if None, uses current figure)
    dpi (int, optional): resolution for static images and base resolution for interactive exports (default: 600)
    PDB_pt (str, optional): File path to PDB file for Mol* visualization (Default: None)
    icon (str, optional): html file icon (Default: python logo)
    """
    if file is None or dir is None:
        return

    mkdir(dir)  # Ensure output directory exists

    if fig is None:
        fig = plt.gcf()

    try:
        fig.tight_layout()
    except Exception:
        pass

    ext = file.split('.')[-1].lower()
    if ext not in ('html', 'json'):  # Save static image
        plt.savefig(
            fname=os.path.join(dir, file),
            dpi=dpi if dpi > 0 else 600,
            bbox_inches='tight',
            format=ext,
        )
    else:
        # Interactive exports via mpld3
        original_dpi = fig.get_dpi()
        fig.set_dpi(dpi if dpi > 0 else 150)

        if ext == 'html':
            if PDB_pt is None:
                # Plain mpld3 HTML
                mpld3.save_html(fig, os.path.join(dir, file))
            else:
                # Copy molstar JS and CSS paths
                sub_dir = os.path.join(dir,file[:-5])
                mkdir(sub_dir)
                with pkg_resources.path(molstar_pkg, "molstar.js") as js_path:
                    shutil.copy(js_path, Path(sub_dir) / "molstar.js")
                with pkg_resources.path(molstar_pkg, "molstar.css") as css_path:
                    shutil.copy(css_path, Path(sub_dir) / "molstar.css")
                with pkg_resources.path(icon_pkg, f"{icon}.svg") as svg_path:
                    shutil.copy(svg_path, Path(sub_dir) / f"{icon}.svg")
                
                # Build a combined HTML: mpld3 plot (left) + Mol* viewer (right)
                if len(PDB_pt) == 4:
                    for PDB_file in os.listdir(os.path.expanduser('~/.config/edms/PDB/')):
                        if PDB_pt in PDB_file:
                            PDB_pt = f'{os.path.expanduser("~/.config/edms/PDB")}/{PDB_file}'
                            break

                export_mpld3_molstar_html(
                    fig=fig,
                    dir=dir,
                    file=file,
                    pdb_url=PDB_pt,
                    icon=icon
                )
                
        elif ext == 'json':
            fig_dict = mpld3.fig_to_dict(fig)
            with open(os.path.join(dir, file), 'w') as f:
                json.dump(fig_dict, f)

        fig.set_dpi(original_dpi)

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

def formatter(typ: str, ax, df: pd.DataFrame, x: str, y: str, cols: str, file: str, dir: str,
              title: str, title_size: int, title_weight: str, title_font: str,
              x_axis: str, x_axis_size: int, x_axis_weight: str, x_axis_font: str, x_axis_scale: str, x_axis_dims: tuple, x_axis_pad: int, x_ticks_size: int, x_ticks_rot: int, x_ticks_font: str, x_ticks: list,
              y_axis: str, y_axis_size: int, y_axis_weight: str, y_axis_font: str, y_axis_scale: str, y_axis_dims: tuple, y_axis_pad: int, y_ticks_size: int, y_ticks_rot: int, y_ticks_font: str, y_ticks: list,
              legend_title: str, legend_title_size: int, legend_size: int, legend_bbox_to_anchor: tuple, legend_loc: str, legend_items: tuple, legend_ncol: int,
              legend_columnspacing: int=-4, legend_handletextpad: float=0.5, legend_labelspacing: float=0.5, legend_borderpad: float=0.5, legend_handlelength: float=0.5, 
              dpi: int = 0, show: bool = True, space_capitalize: bool = True, PDB_pt: str = None, icon: str = 'python', cats_ord: list = None) -> None:
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
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_scale (str, optional): x-axis scale linear, log, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_axis_pad (int, optional): x-axis padding
    x_ticks_size (int, optional): x-axis ticks font size
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-axis ticks font family
    x_ticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_axis_pad (int, optional): y-axis padding
    y_ticks_size (int, optional): y-axis ticks font size
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks_font (str, optional): y-axis ticks font family
    y_ticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_items (tuple, optional): legend items range (start, end)
    legend_ncol (tuple, optional): # of columns
    legend_columnspacing (int, optional): spacing between columns (Default: -4)
    legend_handletextpad (float, optional): spacing between handle and text (Default: 0.5)
    legend_labelspacing (float, optional): vertical spacing between labels (Default: 0.5)
    legend_borderpad (float, optional): border padding (Default: 0.5)
    legend_handlelength (float, optional): length of legend handle (Default: 0.5)
    html (bool, optional): if True, also save an interactive HTML version of the figure using mpld3 (Default: False)
    dpi (int, optional): figure dpi (Default: 600 for non-HTML, 150 for HTML)
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    PDB_pt (str, optional): Path to PDB file for mol* visualization; only intended for fastq plots (Default: None)
    icon (str, optional): html file icon (Default: python logo)
    cats_ord (list, optional): order of categories for categorical plots

    Dependencies: os, matplotlib, seaborn, io, re_un_cap(), & round_up_pow_10()
    '''
    # Define plot types
    scats = ['scat', 'line', 'line_scat']
    cats = ['bar', 'box', 'violin', 'swarm', 'strip', 'point', 'count', 'bar_strip', 'box_strip', 'violin_strip','bar_swarm', 'box_swarm', 'violin_swarm']
    dists = ['hist', 'kde', 'hist_kde','rid']
    heats = ['ht']

    # Determine if file is html
    if file is not None:
        is_html = file.endswith('.html')
    else:
        is_html = False
        
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
        plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight,fontfamily=x_axis_font, labelpad=x_axis_pad)
        if x!='':
            if df[x].apply(lambda row: isinstance(row, (int, float))).all()==True: # Check that x column is numeric
                plt.xscale(x_axis_scale)
                if (x_axis_dims==(0,0))&(x_axis_scale=='log'): plt.xlim(round_down_pow_10(min(df[x])),round_up_pow_10(max(df[x])))
                elif x_axis_dims==(0,0): print('Default x axis dimensions.')
                else: plt.xlim(x_axis_dims[0],x_axis_dims[1])
            else:
                if is_html: # put one tick per category, centered on each bar row
                    ls = df[x].unique().tolist() if cats_ord is None else cats_ord
                    ax.set_xticks(np.arange(len(ls)))
                    ax.set_xticklabels(ls,
                                    rotation=x_ticks_rot,
                                    fontfamily=x_ticks_font,
                                    fontsize=x_ticks_size)
                    ax.set_xlim(-0.5, len(ls) - 0.5)  # avoid extra half-step ticks'''
        if x_ticks==[]: 
            if x_ticks_rot==0: plt.xticks(rotation=x_ticks_rot,ha='center',va='top',fontfamily=x_ticks_font,fontsize=x_ticks_size)
            elif x_ticks_rot == 90: plt.xticks(rotation=x_ticks_rot,ha='right',va='center',fontfamily=x_ticks_font,fontsize=x_ticks_size)
            else: plt.xticks(rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font,fontsize=x_ticks_size)
        else: 
            if x_ticks_rot==0: plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot, ha='center',va='top',fontfamily=x_ticks_font,fontsize=x_ticks_size)
            elif x_ticks_rot == 90: plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot,ha='right',va='center',fontfamily=x_ticks_font,fontsize=x_ticks_size)
            else: plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font,fontsize=x_ticks_size)

        # Set y axis
        if y_axis=='': 
            if space_capitalize: y_axis=re_un_cap(y)
            else: y_axis=y
        plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font, labelpad=y_axis_pad)
        if y!='':
            if df[y].apply(lambda row: isinstance(row, (int, float))).all()==True: # Check that y column is numeric
                plt.yscale(y_axis_scale)
                if (y_axis_dims==(0,0))&(y_axis_scale=='log'): plt.ylim(round_down_pow_10(min(df[y])),round_up_pow_10(max(df[y])))
                elif y_axis_dims==(0,0): print('Default y axis dimensions.')
                else: plt.ylim(y_axis_dims[0],y_axis_dims[1])
            else:
                if is_html: # put one tick per category, centered on each bar row
                    ls = df[y].unique().tolist() if cats_ord is None else cats_ord
                    ax.set_yticks(np.arange(len(ls)))
                    ax.set_yticklabels(ls,
                                    rotation=y_ticks_rot,
                                    fontfamily=y_ticks_font,
                                    fontsize=y_ticks_size)
                    ax.set_ylim(-0.5, len(ls) - 0.5)  # avoid extra half-step ticks
        if y_ticks==[]: plt.yticks(rotation=y_ticks_rot,fontfamily=y_ticks_font, fontsize=y_ticks_size)
        else: plt.yticks(ticks=y_ticks,labels=y_ticks,rotation=y_ticks_rot,fontfamily=y_ticks_font, fontsize=y_ticks_size)

        # Set legend
        if cols is None: print('No legend because cols was not specified.')
        else:
            if legend_title=='': legend_title=cols
            if legend_items==(0,0) and typ not in dists:
                if is_html:
                    if legend_bbox_to_anchor == (1,1): legend_bbox_to_anchor = (-0.1,-0.15)
                    ax.legend(title=legend_title, 
                                title_fontsize=legend_title_size, fontsize=legend_size,
                                bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol,
                                columnspacing=legend_columnspacing,    # space between columns
                                handletextpad=legend_handletextpad,    # space between marker and text
                                labelspacing=legend_labelspacing,      # vertical space between entries
                                borderpad=legend_borderpad,            # padding inside legend box
                                handlelength=legend_handlelength)      # marker length
                else:
                    ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                            bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol) # Move legend to the right of the graph
            elif typ not in dists:
                handles, labels = ax.get_legend_handles_labels()
                if is_html:
                    if legend_bbox_to_anchor == (1,1): legend_bbox_to_anchor = (-0.1,-0.15)
                    ax.legend(title=legend_title, 
                                title_fontsize=legend_title_size, fontsize=legend_size,
                                bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol,
                                handles=handles[legend_items[0]:legend_items[1]],labels=labels[legend_items[0]:legend_items[1]], # Only retains specified labels
                                columnspacing=legend_columnspacing,    # space between columns
                                handletextpad=legend_handletextpad,    # space between marker and text
                                labelspacing=legend_labelspacing,      # vertical space between entries
                                borderpad=legend_borderpad,            # padding inside legend box
                                handlelength=legend_handlelength)      # marker length
                else:
                    ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                            bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol, # Move right of the graph
                            handles=handles[legend_items[0]:legend_items[1]],labels=labels[legend_items[0]:legend_items[1]]) # Only retains specified labels
            else: move_dist_legend(ax,legend_loc,legend_title_size,legend_size,legend_bbox_to_anchor,legend_ncol)

    # Save & show fig
    save_fig(file=file, dir=dir, fig=ax.figure, dpi=dpi, PDB_pt=PDB_pt, icon=icon)
    if show:
        ext = file.split('.')[-1].lower() if file is not None else ''
        if ext not in ('html', 'json'):
            plt.show()
        else: 
            mpld3.show()

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
    elif palette_or_cmap in plt.colormaps() or isinstance(palette_or_cmap, mcolors.Colormap): # Check if cmap is a valid matplotlib color map or custom colormap
        cmap = cm.get_cmap(palette_or_cmap) # Get the color map
        return mcolors.ListedColormap([cmap(i) for i in range(cmap.N)] * repeats) # Breaks the color map into a repeated list
    else:
        print(f'{cmap} is not a valid matplotlib color map and did not apply repeat.')
        return cmap

# Graph methods
def scat(typ: str, df: pd.DataFrame | str, x: str, y: str, cols: str = None, cols_ord: list = None, cols_exclude: list | str = None, 
         stys: str = None, stys_order: list = [], mark_order: list = [], label: str | None = None,
         file: str = None, dir: str = None, palette_or_cmap: str = 'colorblind', alpha: float = 1.0, edgecol: str = 'black',
         figsize: tuple = (5, 5), title: str = '', title_size: int = 18, title_weight: str = 'bold', title_font: str = 'Arial',
         x_axis: str = '', x_axis_size: int = 12, x_axis_weight: str = 'bold', x_axis_font: str = 'Arial', x_axis_scale: str = 'linear', x_axis_dims: tuple = (0, 0), x_axis_pad: int = None, x_ticks_size: int = 9, x_ticks_rot: int = 0, x_ticks_font: str = 'Arial', x_ticks: list = [],
         y_axis: str = '', y_axis_size: int = 12, y_axis_weight: str = 'bold', y_axis_font: str = 'Arial', y_axis_scale: str = 'linear', y_axis_dims: tuple = (0, 0), y_axis_pad: int = None, y_ticks_size: int = 9, y_ticks_rot: int = 0, y_ticks_font: str = 'Arial', y_ticks: list = [],
         legend_title: str = '', legend_title_size: int = 12, legend_size: int = 9, legend_bbox_to_anchor: tuple = (1, 1), legend_loc: str = 'upper left', legend_items: tuple = (0, 0), legend_ncol: int = 1,
         legend_columnspacing: int=0, legend_handletextpad: float=0.5, legend_labelspacing: float=0.5, legend_borderpad: float=0.5, legend_handlelength: float=1, legend_size_html_multiplier: float=1.0,
         dpi: int = 0, show: bool = True, space_capitalize: bool = True, **kwargs):
    ''' 
    scat(): creates scatter plot related graphs

    Parameters:
    typ (str): plot type (scat, line, line_scat)
    df (dataframe | str): pandas dataframe (or file path)
    x (str): x-axis column name
    y (str): y-axis column name
    cols (str, optional): color column name
    cols_ord (list, optional): color column values order
    cols_exclude (list | str, optional): color column values exclude
    stys (str, optional): styles column name
    stys_order (list, optional): styles order
    mark_order (list, optional): marker order
    label (str, optional): column name for point labels; static text for images, interactive tooltips for HTML
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    alpha (float, optional): Alpha (transparency) for scatter points (0 to 1)
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
    x_axis_pad (int, optional): x-axis padding
    x_ticks_size (int, optional): x-axis tick font size
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-ticks font
    x_ticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_axis_pad (int, optional): y-axis padding
    y_ticks_size (int, optional): y-axis tick font size
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks_font (str, optional): y_ticks font
    y_ticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    legend_columnspacing (int, optional): space between columns (Default: 0; only for html plots)
    legend_handletextpad (float, optional): space between marker and text (Default: 0.5; only for html plots)
    legend_labelspacing (float, optional): vertical space between entries (Default: 0.5; only for html plots)
    legend_borderpad (float, optional): padding inside legend box (Default: 0.5; only for html plots)
    legend_handlelength (float, optional): marker length (Default: 1; only for html plots)
    legend_size_html_multiplier (float, optional): multiplier for legend font size in html plots (Default: 1.0)
    dpi (int, optional): figure dpi (Default: 600 for non-HTML, 150 for HTML)
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
    elif palette_or_cmap in plt.colormaps() or isinstance(palette_or_cmap, mcolors.Colormap): 
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
        if typ=='scat': 
            sns.scatterplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, style_order=stys_order, markers=mark_order, edgecolor=edgecol, alpha=alpha, palette=palette, ax=ax, **kwargs)
        elif typ=='line': 
            sns.lineplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, style_order=stys_order, markers=mark_order, palette=palette, ax=ax, **kwargs)
        elif typ=='line_scat':
            sns.lineplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, style_order=stys_order, markers=mark_order, palette=palette, ax=ax, **kwargs)  
            sns.scatterplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, style_order=stys_order, markers=mark_order, edgecolor=edgecol, alpha=alpha, palette=palette, ax=ax, **kwargs)
        else:
            print("Invalid type! scat, line, or line_scat")
            return
    elif cols is not None:
        if typ=='scat': 
            sns.scatterplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, alpha=alpha, palette=palette, ax=ax, **kwargs)
        elif typ=='line': 
            sns.lineplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, ax=ax, palette=palette, **kwargs)
        elif typ=='line_scat':
            sns.lineplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, palette=palette, ax=ax, **kwargs)  
            sns.scatterplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, alpha=alpha, palette=palette, ax=ax, **kwargs)
        else:
            print("Invalid type! scat, line, or line_scat")
            return
    elif stys is not None:
        if typ=='scat': 
            sns.scatterplot(data=df, x=x, y=y, style=stys, style_order=stys_order, markers=mark_order,  edgecolor=edgecol, palette=palette, alpha=alpha, ax=ax, **kwargs)
        elif typ=='line': 
            sns.lineplot(data=df, x=x, y=y, style=stys, style_order=stys_order, markers=mark_order, palette=palette, ax=ax, **kwargs)
        elif typ=='line_scat':
            sns.lineplot(data=df, x=x, y=y, style=stys, style_order=stys_order, markers=mark_order, palette=palette, ax=ax, **kwargs)  
            sns.scatterplot(data=df, x=x, y=y, style=stys, style_order=stys_order, markers=mark_order, edgecolor=edgecol, palette=palette, alpha=alpha, ax=ax, **kwargs)
        else:
            print("Invalid type! scat, line, or line_scat")
            return
    else:
        if typ=='scat': 
            sns.scatterplot(data=df, x=x, y=y, edgecolor=edgecol, palette=palette, alpha=alpha, ax=ax, **kwargs)
        elif typ=='line': 
            sns.lineplot(data=df, x=x, y=y, palette=palette, ax=ax, **kwargs)
        elif typ=='line_scat':
            sns.lineplot(data=df, x=x, y=y, palette=palette, ax=ax, **kwargs)  
            sns.scatterplot(data=df, x=x, y=y, edgecolor=edgecol, palette=palette, alpha=alpha, ax=ax, **kwargs)
        else:
            print("Invalid type! scat, line, or line_scat")
            return
    
    # Determine if file is html
    if file is not None:
        is_html = file.endswith('.html')
    else:
        is_html = False
    
    # Match title fontsize for html plots
    if is_html:
        x_axis_size=title_size
        y_axis_size=title_size
        x_ticks_size=title_size
        y_ticks_size=title_size
        legend_title_size=title_size
        legend_size=title_size*legend_size_html_multiplier

    # Optional point labels: static for images, interactive for HTML
    if label is not None and label in df.columns:
        if is_html:
            # Use invisible points to carry click-to-toggle tooltips in HTML
            pts = ax.scatter(
                df[x],
                df[y],
                s=20,
                alpha=0
            )
            labels_list = df[label].fillna("").astype(str).tolist()
            tooltip = SafeHTMLTooltip(pts, labels_list)
            clicker = ClickTooltip(pts, labels_list)
            mpld3.plugins.connect(fig, tooltip, clicker)

        else:
            # For static images, draw permanent text labels
            for i, txt in enumerate(df[label]):
                ax.text(
                    df.iloc[i][x],
                    df.iloc[i][y],
                    txt
                )

    formatter(typ=typ, ax=ax, df=df, x=x, y=y, cols=cols, file=file, dir=dir, 
              title=title, title_size=title_size, title_weight=title_weight, title_font=title_font,
              x_axis=x_axis, x_axis_size=x_axis_size, x_axis_weight=x_axis_weight, x_axis_font=x_axis_font, x_axis_scale=x_axis_scale, x_axis_dims=x_axis_dims, x_axis_pad=x_axis_pad, x_ticks_size=x_ticks_size, x_ticks_rot=x_ticks_rot, x_ticks_font=x_ticks_font, x_ticks=x_ticks,
              y_axis=y_axis, y_axis_size=y_axis_size, y_axis_weight=y_axis_weight, y_axis_font=y_axis_font, y_axis_scale=y_axis_scale, y_axis_dims=y_axis_dims, y_axis_pad=y_axis_pad, y_ticks_size=y_ticks_size, y_ticks_rot=y_ticks_rot, y_ticks_font=y_ticks_font, y_ticks=y_ticks,
              legend_title=legend_title, legend_title_size=legend_title_size, legend_size=legend_size, legend_bbox_to_anchor=legend_bbox_to_anchor, legend_loc=legend_loc, legend_items=legend_items, legend_ncol=legend_ncol,
              legend_columnspacing=legend_columnspacing, legend_handletextpad=legend_handletextpad, legend_labelspacing=legend_labelspacing, legend_borderpad=legend_borderpad, legend_handlelength=legend_handlelength,
              dpi=dpi, show=show, space_capitalize=space_capitalize, icon='scatter')

def cat(typ: str, df: pd.DataFrame | str, x: str = '', y: str = '', cats_ord: list = None, cats_exclude: list|str = None, cols: str = None, cols_ord: list = None, cols_exclude: list | str = None, line: float = None,
        file: str = None, dir: str = None, palette_or_cmap: str = 'colorblind', alpha: float = 1.0, dodge: bool = False, jitter: bool = True, size: float = 5, edgecol: str = 'black', lw: int = 1, errorbar: str = 'sd', errwid: int = 1, errcap: float = 0.1,
        figsize: tuple = (5, 5), title: str = '', title_size: int = 18, title_weight: str = 'bold', title_font: str = 'Arial',
        x_axis: str = '', x_axis_size: int = 12, x_axis_weight: str = 'bold', x_axis_font: str = 'Arial', x_axis_scale: str = 'linear', x_axis_dims: tuple = (0, 0), x_axis_pad: int = None, x_ticks_size: int = 9, x_ticks_rot: int = 0, x_ticks_font: str = 'Arial', x_ticks: list = [],
        y_axis: str = '', y_axis_size: int = 12, y_axis_weight: str = 'bold', y_axis_font: str = 'Arial', y_axis_scale: str = 'linear', y_axis_dims: tuple = (0, 0), y_axis_pad: int = None, y_ticks_size: int = 9, y_ticks_rot: int = 0, y_ticks_font: str = 'Arial', y_ticks: list = [],
        legend_title: str = '', legend_title_size: int = 12, legend_size: int = 9, legend_bbox_to_anchor: tuple = (1, 1), legend_loc: str = 'upper left', legend_items: tuple = (0, 0), legend_ncol: int = 1,
        legend_columnspacing: int=0, legend_handletextpad: float=0.5, legend_labelspacing: float=0.5, legend_borderpad: float=0.5, legend_handlelength: float=1, legend_size_html_multiplier: float=1.0,
        dpi: int = 0, show: bool = True, space_capitalize: bool = True, PDB_pt: str = None, **kwargs):
    ''' 
    cat(): creates categorical graphs

    Parameters:
    typ (str): plot type (bar, box, violin, swarm, strip, point, count, bar_strip, box_strip, violin_strip, bar_swarm, box_swarm, violin_swarm)
    df (dataframe | str): pandas dataframe (or file path)
    x (str, optional): x-axis column name
    y (str, optional): y-axis column name
    cats_ord (list, optional): category column values order (x- or y-axis)
    cats_exclude (list | str, optional): category column values exclude (x- or y-axis)
    cols (str, optional): color column name
    cols_ord (list, optional): color column values order
    cols_exclude (list | str, optional): color column values exclude
    line (float, optional): add horizontal line at y value or vertical line at x value
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    alpha (float, optional): Alpha (transparency) for scatter points (0 to 1)
    dodge (bool, optional): whether to separate points by color category
    jitter (bool, optional): whether to add jitter to points (for strip plots)
    size (float, optional): Marker size for scatter points
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
    x_axis_pad (int, optional): x-axis padding
    x_ticks_size (int, optional): x-axis ticks font size
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-ticks font
    x_ticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_axis_pad (int, optional): y-axis padding
    y_ticks_size (int, optional): y-axis ticks font size
    y_ticks_font (str, optional): y_ticks font
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    legends_items (tuple, optional): legend items to show (start, end)
    legend_columnspacing (int, optional): space between columns (Default: 0; only for html plots)
    legend_handletextpad (float, optional): space between marker and text (Default: 0.5; only for html plots)
    legend_labelspacing (float, optional): vertical space between entries (Default: 0.5; only for html plots)
    legend_borderpad (float, optional): padding inside legend box (Default: 0.5; only for html plots)
    legend_handlelength (float, optional): marker length (Default: 1; only for html plots)
    legend_size_html_multiplier (float, optional): multiplier for legend font size in html plots (Default: 1.0)
    dpi (int, optional): figure dpi (Default: 600 for non-HTML, 150 for HTML)
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    PDB_pt (str, optional): Path to PDB file for mol* visualization; only intended for fastq plots (Default: None)
    
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), & round_up_pow_10()
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)
    
    # Match title fontsize for html plots
    if file is not None:
        if file.endswith('.html')==True:
            x_axis_size=title_size
            y_axis_size=title_size
            x_ticks_size=title_size
            y_ticks_size=title_size
            legend_title_size=title_size
            legend_size=title_size*legend_size_html_multiplier

    # Omit excluded data
    if type(cats_exclude)==list: # Category exclusion
        for exclude in cats_exclude: 
            if x!='': df=df[df[x]!=exclude]
            if y!='': df=df[df[y]!=exclude]
    elif type(cats_exclude)==str: 
        if x!='': df=df[df[x]!=cats_exclude]
        if y!='': df=df[df[y]!=cats_exclude]

    if type(cols_exclude)==list: # Color exclusion
        for exclude in cols_exclude: df=df[df[cols]!=exclude]
    elif type(cols_exclude)==str: df=df[df[cols]!=cols_exclude]

    # Set color scheme and category column order
    color_palettes = ["deep", "muted", "bright", "pastel", "dark", "colorblind", "husl", "hsv", "Paired", "Set1", "Set2", "Set3", "tab10", "tab20"] # List of common Seaborn palettes
    if palette_or_cmap in color_palettes: palette = palette_or_cmap
    elif palette_or_cmap in plt.colormaps() or isinstance(palette_or_cmap, mcolors.Colormap): 
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

        if typ=='bar': 
            sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, err_kws={'color':edgecol, 'linewidth':errwid}, capsize=errcap, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='box': 
            sns.boxplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='violin': 
            sns.violinplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='swarm': 
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='strip': 
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, jitter=jitter, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='point': 
            sns.pointplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, err_kws={'linewidth':errwid}, capsize=errcap, hue=cols, hue_order=cols_ord, palette=palette, ax=ax, **kwargs)
        elif typ=='count': 
            if (x!='')&(y!=''):
                print('Cannot make countplot with both x and y specified.')
                return
            elif x!='': 
                sns.countplot(data=df, x=x, order=cats_ord, hue=cols, hue_order=cols_ord, palette=palette, ax=ax, **kwargs)
            elif y!='': 
                sns.countplot(data=df, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, palette=palette, ax=ax, **kwargs)
            else:
                print('Cannot make countplot without x or y specified.')
                return
        elif typ=='bar_strip':
            sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, err_kws={'color':edgecol, 'linewidth':errwid}, capsize=errcap, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, jitter=jitter, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='box_strip':
            sns.boxplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, jitter=jitter, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_strip':
            sns.violinplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, jitter=jitter, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='bar_swarm':
            sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, err_kws={'color':edgecol, 'linewidth':errwid}, capsize=errcap, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='box_swarm':
            sns.boxplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_swarm':
            sns.violinplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, size=size, palette=palette, ax=ax, **kwargs)
        else:
            print('Invalid type! bar, box, violin, swarm, strip, point, count, bar_strip, box_strip, violin_strip, bar_swarm, box_swarm, violin_swarm')
            return

    else: # Cols was not specified
        
        if typ=='bar': 
            sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, err_kws={'color':edgecol, 'linewidth':errwid}, capsize=errcap, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='box': 
            sns.boxplot(data=df, x=x, y=y, order=cats_ord, linewidth=lw, ax=ax, palette=palette, **kwargs)
        elif typ=='violin': 
            sns.violinplot(data=df, x=x, y=y, order=cats_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='swarm': 
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='strip': 
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, jitter=jitter, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='point': 
            sns.pointplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, err_kws={'linewidth':errwid}, capsize=errcap, palette=palette, ax=ax, **kwargs)
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
            sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, err_kws={'color':edgecol, 'linewidth':errwid}, capsize=errcap, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, jitter=jitter, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='box_strip':
            sns.boxplot(data=df, x=x, y=y, order=cats_ord, linewidth=lw, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, jitter=jitter, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_strip':
            sns.violinplot(data=df, x=x, y=y, order=cats_ord, edgecolor=edgecol, linewidth=lw, ax=ax, palette=palette, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, jitter=jitter, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='bar_swarm':
            sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, err_kws={'color':edgecol, 'linewidth':errwid}, capsize=errcap, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='box_swarm':
            sns.boxplot(data=df, x=x, y=y, order=cats_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, size=size, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_swarm':
            sns.violinplot(data=df, x=x, y=y, order=cats_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, alpha=alpha, linewidth=lw, dodge=dodge, size=size, palette=palette, ax=ax, **kwargs)
        else:
            print('Invalid type! bar, box, violin, swarm, strip, point, count, bar_strip, box_strip, violin_strip, bar_swarm, box_swarm, violin_swarm')
            return
    
    # Add line if specified
    if line is not None:
        if df[x].apply(lambda row: isinstance(row, str)).all()==True: # x column is categorical; y is numeric
            ax.axhline(y=line, color='black', linestyle='--', linewidth=lw)
        elif df[y].apply(lambda row: isinstance(row, str)).all()==True: # y column is categorical; x is numeric
            ax.axvline(x=line, color='black', linestyle='--', linewidth=lw)
        else:
            print('Cannot add line: both x and y are numeric.')


    formatter(typ=typ, ax=ax, df=df, x=x, y=y, cols=cols, file=file, dir=dir, 
              title=title, title_size=title_size, title_weight=title_weight, title_font=title_font,
              x_axis=x_axis, x_axis_size=x_axis_size, x_axis_weight=x_axis_weight, x_axis_font=x_axis_font, x_axis_scale=x_axis_scale, x_axis_dims=x_axis_dims, x_axis_pad=x_axis_pad, x_ticks_size=x_ticks_size, x_ticks_rot=x_ticks_rot, x_ticks_font=x_ticks_font, x_ticks=x_ticks,
              y_axis=y_axis, y_axis_size=y_axis_size, y_axis_weight=y_axis_weight, y_axis_font=y_axis_font, y_axis_scale=y_axis_scale, y_axis_dims=y_axis_dims, y_axis_pad=y_axis_pad, y_ticks_size=y_ticks_size, y_ticks_rot=y_ticks_rot, y_ticks_font=y_ticks_font, y_ticks=y_ticks,
              legend_title=legend_title, legend_title_size=legend_title_size, legend_size=legend_size, legend_bbox_to_anchor=legend_bbox_to_anchor, legend_loc=legend_loc, legend_items=legend_items, legend_ncol=legend_ncol,
              legend_columnspacing=legend_columnspacing, legend_handletextpad=legend_handletextpad, legend_labelspacing=legend_labelspacing, legend_borderpad=legend_borderpad, legend_handlelength=legend_handlelength,
              dpi=dpi, show=show, space_capitalize=space_capitalize, PDB_pt=PDB_pt, icon='cat', cats_ord=cats_ord)

def dist(typ: str, df: pd.DataFrame | str, x: str, cols: str = None, cols_ord: list = None, cols_exclude: list | str = None, bins: int = 40, log10_low: int = 0,
         file: str = None, dir: str = None, palette_or_cmap: str = 'colorblind', edgecol: str = 'black', lw: int = 1, ht: float = 1.5, asp: int = 5, tp: float = .8, hs: int = 0, despine: bool = False,
         figsize: tuple = (5, 5), title: str = '', title_size: int = 18, title_weight: str = 'bold', title_font: str = 'Arial',
         x_axis: str = '', x_axis_size: int = 12, x_axis_weight: str = 'bold', x_axis_font: str = 'Arial', x_axis_scale: str = 'linear', x_axis_dims: tuple = (0, 0), x_axis_pad: int = None, x_ticks_size: int = 9, x_ticks_rot: int = 0, x_ticks_font: str = 'Arial', x_ticks: list = [],
         y_axis: str = '', y_axis_size: int = 12, y_axis_weight: str = 'bold', y_axis_font: str = 'Arial', y_axis_scale: str = 'linear', y_axis_dims: tuple = (0, 0), y_axis_pad: int = None, y_ticks_size: int = 9, y_ticks_rot: int = 0, y_ticks_font: str = 'Arial', y_ticks: list = [],
         legend_title: str = '', legend_title_size: int = 12, legend_size: int = 9, legend_bbox_to_anchor: tuple = (1, 1), legend_loc: str = 'upper left', legend_items: tuple = (0, 0), legend_ncol: int = 1,
         dpi: int = 0, show: bool = True, space_capitalize: bool = True, **kwargs):
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
    x_axis_pad (int, optional): x-axis label padding
    x_ticks_size (int, optional): x-axis ticks font size
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-ticks font
    x_ticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_axis_pad (int, optional): y-axis label padding
    y_ticks_size (int, optional): y-axis ticks font size
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks_font (str, optional): y_ticks font
    y_ticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    dpi (int, optional): figure dpi (Default: 600 for non-HTML, 150 for HTML)
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
    elif palette_or_cmap in plt.colormaps() or isinstance(palette_or_cmap, mcolors.Colormap): 
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
        formatter(typ=typ, ax=ax, df=df, x=x, y=y, cols=cols, file=file, dir=dir,
                  title=title, title_size=title_size, title_weight=title_weight, title_font=title_font,
                  x_axis=x_axis, x_axis_size=x_axis_size, x_axis_weight=x_axis_weight, x_axis_font=x_axis_font, x_axis_scale=x_axis_scale, x_axis_dims=x_axis_dims, x_axis_pad=x_axis_pad, x_ticks_size=x_ticks_size, x_ticks_rot=x_ticks_rot, x_ticks_font=x_ticks_font, x_ticks=x_ticks,
                  y_axis=y_axis, y_axis_size=y_axis_size, y_axis_weight=y_axis_weight, y_axis_font=y_axis_font, y_axis_scale=y_axis_scale, y_axis_dims=y_axis_dims, y_axis_pad=y_axis_pad, y_ticks_size=y_ticks_size, y_ticks_rot=y_ticks_rot, y_ticks_font=y_ticks_font, y_ticks=y_ticks,
                  legend_title=legend_title, legend_title_size=legend_title_size, legend_size=legend_size, legend_bbox_to_anchor=legend_bbox_to_anchor, legend_loc=legend_loc, legend_items=legend_items, legend_ncol=legend_ncol,
                  dpi=dpi, show=show, space_capitalize=space_capitalize, icon='histogram')
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
        formatter(typ=typ, ax=ax, df=df, x=x, y=y, cols=cols, file=file, dir=dir, 
                  title=title, title_size=title_size, title_weight=title_weight, title_font=title_font,
                  x_axis=x_axis, x_axis_size=x_axis_size, x_axis_weight=x_axis_weight, x_axis_font=x_axis_font, x_axis_scale=x_axis_scale, x_axis_dims=x_axis_dims, x_axis_pad=x_axis_pad, x_ticks_size=x_ticks_size, x_ticks_rot=x_ticks_rot, x_ticks_font=x_ticks_font, x_ticks=x_ticks,
                  y_axis=y_axis, y_axis_size=y_axis_size, y_axis_weight=y_axis_weight, y_axis_font=y_axis_font, y_axis_scale=y_axis_scale, y_axis_dims=y_axis_dims, y_axis_pad=y_axis_pad, y_ticks_size=y_ticks_size, y_ticks_rot=y_ticks_rot, y_ticks_font=y_ticks_font, y_ticks=y_ticks,
                  legend_title=legend_title, legend_title_size=legend_title_size, legend_size=legend_size, legend_bbox_to_anchor=legend_bbox_to_anchor, legend_loc=legend_loc, legend_items=legend_items, legend_ncol=legend_ncol,
                  dpi=dpi, show=show, space_capitalize=space_capitalize, icon='histogram')
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
        formatter(typ=typ, ax=ax, df=df, x=x, y=y, cols=cols, file=file, dir=dir, 
                  title=title, title_size=title_size, title_weight=title_weight, title_font=title_font,
                  x_axis=x_axis, x_axis_size=x_axis_size, x_axis_weight=x_axis_weight, x_axis_font=x_axis_font, x_axis_scale=x_axis_scale, x_axis_dims=x_axis_dims, x_axis_pad=x_axis_pad, x_ticks_size=x_ticks_size, x_ticks_rot=x_ticks_rot, x_ticks_font=x_ticks_font, x_ticks=x_ticks,
                  y_axis=y_axis, y_axis_size=y_axis_size, y_axis_weight=y_axis_weight, y_axis_font=y_axis_font, y_axis_scale=y_axis_scale, y_axis_dims=y_axis_dims, y_axis_pad=y_axis_pad, y_ticks_size=y_ticks_size, y_ticks_rot=y_ticks_rot, y_ticks_font=y_ticks_font, y_ticks=y_ticks,
                  legend_title=legend_title, legend_title_size=legend_title_size, legend_size=legend_size, legend_bbox_to_anchor=legend_bbox_to_anchor, legend_loc=legend_loc, legend_items=legend_items, legend_ncol=legend_ncol,
                  dpi=dpi, show=show, space_capitalize=space_capitalize, icon='histogram')
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
        save_fig(file=file, dir=dir, fig=g.figure, dpi=dpi, icon='histogram')
        if show:
            ext = file.split('.')[-1].lower()  if file is not None else ''
            if ext not in ('html', 'json'):
                plt.show()
            else: 
                mpld3.show()
    else:
        print('Invalid type! hist, kde, hist_kde, rid')
        return

def heat(df: pd.DataFrame | str, x: str = None, y: str = None, vars: str = None, vals: str = None, vals_dims: tuple = None,
         file: str = None, dir: str = None, edgecol: str = 'black', lw: int = 1, annot: bool = False, center: float = None, cmap: str = "Reds", sq: bool = True,
         cbar: bool=True, cbar_label: str=None, cbar_label_size: int=None, cbar_label_weight: str='bold', cbar_tick_size: int=None, cbar_shrink: float=None, cbar_aspect: int=None, cbar_pad: float=None, cbar_orientation: str=None,
         title: str = '', title_size: int = 18, title_weight: str = 'bold', title_font: str = 'Arial', figsize: tuple = (5, 5),
         x_axis: str = '', x_axis_size: int = 12, x_axis_weight: str = 'bold', x_axis_font: str = 'Arial', x_axis_pad: int = None, x_ticks_size: int = 9, x_ticks_rot: int = 45, x_ticks_font: str = 'Arial',
         y_axis: str = '', y_axis_size: int = 12, y_axis_weight: str = 'bold', y_axis_font: str = 'Arial', y_axis_pad: int = None, y_ticks_size: int = 9, y_ticks_rot: int = 0, y_ticks_font: str = 'Arial',
         dpi: int = 0, show: bool = True, space_capitalize: bool = True, **kwargs):
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
    center (float, optional): center value for colormap
    cmap (str, optional): matplotlib color map
    sq (bool, optional): square dimensions (Default: True)
    cbar (bool, optional): show colorbar (Default: True)
    cbar_label (str, optional): colorbar label
    cbar_label_size (int, optional): colorbar label font size
    cbar_label_weight (str, optional): colorbar label bold, normal, & heavy
    cbar_tick_size (int, optional): colorbar tick font size
    cbar_shrink (float, optional): colorbar shrink factor
    cbar_aspect (int, optional): colorbar aspect ratio
    cbar_pad (float, optional): colorbar padding
    cbar_orientation (str, optional): colorbar orientation ('vertical' | 'horizontal')
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    title_font (str, optional): plot title font
    figsize (tuple, optional): figure size per subplot
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_font (str, optional): x-axis font
    x_axis_pad (int, optional): x-axis label padding
    x_ticks_size (int, optional): x-axis ticks font size
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-ticks font
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_axis_pad (int, optional): y-axis label padding
    y_ticks_size (int, optional): y-axis ticks font size
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks_font (str, optional): y-ticks font
    dpi (int, optional): figure dpi (Default: 600 for non-HTML, 150 for HTML)
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), & round_up_pow_10()
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)
    
    # cbar kwargs
    cbar_kws = dict()
    if cbar_label is not None: cbar_kws['label'] = cbar_label
    if cbar_shrink is not None: cbar_kws['shrink'] = cbar_shrink
    if cbar_aspect is not None: cbar_kws['aspect'] = cbar_aspect
    if cbar_pad is not None: cbar_kws['pad'] = cbar_pad
    if cbar_orientation is not None: cbar_kws['orientation'] = cbar_orientation

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
        sns.heatmap(dc[key],annot=annot,cmap=cmap,ax=ax,linecolor=edgecol,linewidths=lw,cbar=cbar,square=sq,vmin=vmin,vmax=vmax,cbar_kws=cbar_kws, **kwargs)
        
        # Title
        if len(list(dc.keys()))>1: ax.set_title(key,fontsize=title_size,fontweight=title_weight,fontfamily=title_font)  # Add title to subplot
        else: ax.set_title(title,fontsize=title_size,fontweight=title_weight,fontfamily=title_font)
        
        # X-axis
        if x_axis=='': 
            if space_capitalize: ax.set_xlabel(re_un_cap(x),fontsize=x_axis_size,fontweight=x_axis_weight,fontfamily=x_axis_font,labelpad=x_axis_pad) # Add x axis label
            else: ax.set_xlabel(x,fontsize=x_axis_size,fontweight=x_axis_weight,fontfamily=x_axis_font,labelpad=x_axis_pad) # Add x axis label
        else: ax.set_xlabel(x_axis,fontsize=x_axis_size,fontweight=x_axis_weight,fontfamily=x_axis_font,labelpad=x_axis_pad)
        
        # Y-axis
        if y_axis=='': 
            if space_capitalize: ax.set_ylabel(re_un_cap(y),fontsize=y_axis_size,fontweight=y_axis_weight,fontfamily=y_axis_font,labelpad=y_axis_pad) # Add y axis label
            else: ax.set_ylabel(y,fontsize=y_axis_size,fontweight=y_axis_weight,fontfamily=y_axis_font,labelpad=y_axis_pad) # Add y axis label
        else: ax.set_ylabel(y_axis,fontsize=y_axis_size,fontweight=y_axis_weight,fontfamily=y_axis_font,labelpad=y_axis_pad)
        
        # Format x ticks
        if x_ticks_rot==0: plt.setp(ax.get_xticklabels(), rotation=x_ticks_rot, ha="center", va='top', rotation_mode="anchor",fontname=x_ticks_font,fontsize=x_ticks_size) 
        elif x_ticks_rot==90: plt.setp(ax.get_xticklabels(), rotation=x_ticks_rot, ha="right", va='center', rotation_mode="anchor",fontname=x_ticks_font,fontsize=x_ticks_size) 
        else: plt.setp(ax.get_xticklabels(), rotation=x_ticks_rot, ha="right", va='top', rotation_mode="anchor",fontname=x_ticks_font,fontsize=x_ticks_size) 
        
        # Format y ticks
        plt.setp(ax.get_yticklabels(), rotation=y_ticks_rot, va='center', ha="right",rotation_mode="anchor",fontname=y_ticks_font,fontsize=y_ticks_size)
        
        # Format cbar
        cbar = ax.collections[0].colorbar
        vmin, vmax = cbar.vmin, cbar.vmax
        if center is None: center = (vmin + vmax) / 2
        cbar.set_label(cbar_label, fontsize=cbar_label_size, fontweight=cbar_label_weight)
        cbar.set_ticks([vmin, center, vmax])
        cbar.ax.tick_params(labelsize=cbar_tick_size)

        # Set white background
        ax.set_facecolor('white')
    
    # Save & show fig
    save_fig(file=file, dir=dir, fig=fig, dpi=dpi, icon='heat')
    if show:
        ext = file.split('.')[-1].lower()  if file is not None else ''
        if ext not in ('html', 'json'):
            plt.show()
        else: 
            mpld3.show()

def stack(df: pd.DataFrame | str, x: str, y: str, cols: str, cutoff_group: str = '', cutoff_value: float = 0, cutoff_keep: bool = True, cols_ord: list = [], x_ord: list = [],
          file: str = None, dir: str = None, palette_or_cmap: str = 'tab20', repeats: int = 1, errcap: int = 4, vertical: bool = True,
          figsize: tuple = (5, 5), title: str = '', title_size: int = 18, title_weight: str = 'bold', title_font: str = 'Arial',
          x_axis: str = '', x_axis_size: int = 12, x_axis_weight: str = 'bold', x_axis_font: str = 'Arial', x_axis_pad: int = None, x_ticks_size: int = 9, x_ticks_rot: int = 0, x_ticks_font: str = 'Arial',
          y_axis: str = '', y_axis_size: int = 12, y_axis_weight: str = 'bold', y_axis_font: str = 'Arial', y_axis_dims: tuple = (0, 0),  y_axis_pad: int = None, y_ticks_size: int = 9, y_ticks_rot: int = 0, y_ticks_font: str = 'Arial',
          legend_title: str = '', legend_title_size: int = 12, legend_size: int = 12, legend_bbox_to_anchor: tuple = (1, 1), legend_loc: str = 'upper left', legend_ncol: int = 1, 
          legend_columnspacing: int=0, legend_handletextpad: float=0.5, legend_labelspacing: float=0.5, legend_borderpad: float=0.5, legend_handlelength: float=1, legend_size_html_multiplier: float=1,
          dpi: int = 0, show: bool = True, space_capitalize: bool = True, PDB_pt: str = None, **kwargs):
    ''' 
    stack(): creates stacked bar plot

    Parameters:
    df (dataframe | str): pandas dataframe (or file path)
    x (str, optional): x-axis column name
    y (str, optional): y-axis column name
    cols (str, optional): color column name
    cutoff_group (str, optional): column name to group by when applying cutoff
    cutoff_value (float, optional): y-axis values needs be greater than (e.g. 0)
    cutoff_keep (bool, optional): keep cutoff group even if below cutoff (Default: True)
    cols_ord (list, optional): color column values order
    x_ord (list, optional): x-axis column values order
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
    x_axis_pad (int, optional): x-axis label padding
    x_ticks_size (int, optional): x-axis ticks font size
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-axis ticks font
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_axis_pad (int, optional): y-axis label padding
    y_ticks_size (int, optional): y-axis ticks font size
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks_font (str, optional): y-ticks font
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    legend_columnspacing (int, optional): spacing between legend columns (Default: 0; only for html plots)
    legend_handletextpad (float, optional): pad between legend handle and text (Default: 0.5; only for html plots)
    legend_labelspacing (float, optional): vertical space between legend entries (Default: 0.5; only for html plots)
    legend_borderpad (float, optional): pad between legend and axes (Default: 0.5; only for html plots)
    legend_handlelength (float, optional): length of legend handle (Default: 1; only for html plots)
    legend_size_html_multiplier (float, optional): legend size multiplier for HTML plots (Default: .75)
    dpi (int, optional): figure dpi (Default: 600 for non-HTML, 150 for HTML)
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    PDB_pt (str, optional): Path to PDB file for mol* visualization; only intended for fastq plots (Default: None)
    
    Dependencies: re, os, pandas, numpy, matplotlib.pyplot, & io
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)
    
    # Determine if HTML output
    if file is not None:
        is_html=file.endswith('.html')
    else:
        is_html=False
    
    # Match title fontsize for html plots
    if is_html:
        if file.endswith('.html')==True:
            x_axis_size=title_size
            y_axis_size=title_size
            x_ticks_size=title_size
            y_ticks_size=title_size
            legend_title_size=title_size
            legend_size=title_size*legend_size_html_multiplier

    # Omit smaller than cutoff and convert it to <cutoff
    if cutoff_group in df.columns and cutoff_value>0: # If cutoff group and value specified
        df_cut=df[df[y]>=cutoff_value]
        if cutoff_keep==True: # Keep cutoff group even if below cutoff
            df_other=df[df[y]<cutoff_value]
            for group in list(df_other[cutoff_group].value_counts().keys()):
                df_temp = df_other[df_other[cutoff_group]==group]
                df_temp[y]=sum(df_temp[y])
                df_temp[cols]=f'<{cutoff_value}'
                df_cut = pd.concat([df_cut,df_temp.iloc[:1]])
    else: # Otherwise use full dataframe
        df_cut=df

    # Make pivot table
    df_pivot=pd.pivot_table(df_cut, index=x, columns=cols, values=y, aggfunc=np.mean)
    df_pivot_err=pd.pivot_table(df_cut, index=x, columns=cols, values=y, aggfunc=np.std)
    if cols_ord!=[]: df_pivot=df_pivot.reindex(columns=cols_ord)
    if x_ord!=[]: df_pivot=df_pivot.reindex(index=x_ord)

    # Make stacked barplot
    if vertical: # orientation
        ax = df_pivot.plot(kind='bar',yerr=df_pivot_err,capsize=errcap, figsize=figsize,colormap=repeat_palette_cmap(palette_or_cmap,repeats),stacked=True,**kwargs)
        fig = ax.get_figure()
        
        # Set x axis
        if x_axis=='': 
            if space_capitalize: x_axis=re_un_cap(x)
            else: x_axis=x
        plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight,fontfamily=x_axis_font,labelpad=x_axis_pad)
        if x_ticks_rot == 0: plt.xticks(rotation=x_ticks_rot,ha='center',va='top',fontfamily=x_ticks_font,fontsize=x_ticks_size)
        elif x_ticks_rot == 90: plt.xticks(rotation=x_ticks_rot,ha='right',va='center',fontfamily=x_ticks_font,fontsize=x_ticks_size)
        else: plt.xticks(rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font,fontsize=x_ticks_size)

        # Set y axis
        if y_axis=='': 
            if space_capitalize: y_axis=re_un_cap(y)
            else: y_axis=y
        plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font,labelpad=y_axis_pad)
        plt.yticks(rotation=y_ticks_rot,fontfamily=y_ticks_font, fontsize=y_ticks_size)

        if y_axis_dims==(0,0): print('Default y axis dimensions.')
        else: plt.ylim(y_axis_dims[0],y_axis_dims[1])

    else: # Horizontal orientation
        ax = df_pivot.plot(kind='barh',xerr=df_pivot_err,capsize=errcap, figsize=figsize,colormap=repeat_palette_cmap(palette_or_cmap,repeats),stacked=True,**kwargs)
        fig = ax.get_figure()

        # Set y axis
        if x_axis=='': 
            if space_capitalize: x_axis=re_un_cap(x)
            else: x_axis=x
        plt.ylabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight,fontfamily=x_axis_font,labelpad=x_axis_pad)
        plt.yticks(rotation=x_ticks_rot,fontfamily=x_ticks_font, fontsize=x_ticks_size)
        
        # Set x axis
        if y_axis=='': 
            if space_capitalize: y_axis=re_un_cap(y)
            else: y_axis=y
        plt.xlabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font,labelpad=y_axis_pad)
        if y_ticks_rot == 0: plt.xticks(rotation=y_ticks_rot,ha='center',va='top',fontfamily=y_ticks_font, fontsize=y_ticks_size)
        elif y_ticks_rot == 90: plt.xticks(rotation=y_ticks_rot,ha='right',va='center',fontfamily=y_ticks_font, fontsize=y_ticks_size)
        else: plt.xticks(rotation=y_ticks_rot,ha='right',fontfamily=y_ticks_font, fontsize=y_ticks_size)

        if y_axis_dims==(0,0): print('Default x axis dimensions.')
        else: plt.xlim(y_axis_dims[0],y_axis_dims[1])
        
    # Set title
    if title=='' and file is not None: title=re_un_cap(".".join(file.split(".")[:-1]))
    plt.title(title, fontsize=title_size, fontweight=title_weight, family=title_font)
    
    # Set legend
    if legend_title=='': 
        if space_capitalize: legend_title=re_un_cap(cols)
        else: legend_title=cols
    if is_html:
        if legend_bbox_to_anchor == (1,1): legend_bbox_to_anchor = (-0.1,-0.15)
        plt.legend(title=legend_title, 
                    title_fontsize=legend_title_size, fontsize=legend_size,
                    bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol,
                    columnspacing=legend_columnspacing,    # space between columns
                    handletextpad=legend_handletextpad,    # space between marker and text
                    labelspacing=legend_labelspacing,      # vertical space between entries
                    borderpad=legend_borderpad,            # padding inside legend box
                    handlelength=legend_handlelength)      # marker length
    else:
        plt.legend(title=legend_title, title_fontsize=legend_title_size, fontsize=legend_size, 
                bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol)

    # Save & show fig; return dataframe
    save_fig(file=file, dir=dir, fig=fig, dpi=dpi, PDB_pt=PDB_pt, icon='bar')
    if show:
        ext = file.split('.')[-1].lower() if file is not None else ''
        if ext not in ('html', 'json'):
            plt.show()
        else:
            mpld3.show(fig)

def vol(df: pd.DataFrame | str, FC: str, pval: str, stys: str = None, size: str | bool = None, size_dims: tuple = None, z_col: str = None, z_var: str = None, label: str = None, stys_order: list = [], mark_order: list = [],
        FC_threshold: float = 1, pval_threshold: float = 1, file: str = None, dir: str = None, color: str = 'lightgray', alpha: float = 0.5, edgecol: str = 'black', vertical: bool = True,
        figsize: tuple = (5, 5), title: str = '', title_size: int = 18, title_weight: str = 'bold', title_font: str = 'Arial',
        x_axis: str = '', x_axis_size: int = 12, x_axis_weight: str = 'bold', x_axis_font: str = 'Arial', x_axis_dims: tuple = (0, 0), x_axis_pad: int = None, x_ticks_size: int = 9, x_ticks_rot: int = 0, x_ticks_font: str = 'Arial', x_ticks: list = [],
        y_axis: str = '', y_axis_size: int = 12, y_axis_weight: str = 'bold', y_axis_font: str = 'Arial', y_axis_dims: tuple = (0, 0), y_axis_pad: int = None, y_ticks_size: int = 9, y_ticks_rot: int = 0, y_ticks_font: str = 'Arial', y_ticks: list = [],
        legend_title: str = '', legend_title_size: int = 12, legend_size: int = 9, legend_bbox_to_anchor: tuple = (1, 1), legend_loc: str = 'upper left', legend_ncol: int = 1, 
        legend_columnspacing: int=-4, legend_handletextpad: float=0.5, legend_labelspacing: float=0.5, legend_borderpad: float=0.5, legend_handlelength: float=0.5, legend_size_html_multiplier: float=0.75,
        display_legend: bool = True, display_labels: str = 'FC & p-value', display_lines: bool = False, display_axis: bool = True, return_df: bool = True, dpi: int = 0, show: bool = True, space_capitalize: bool = True,
        PDB_pt: str = None, # only intended for fastq plots
        **kwargs) -> pd.DataFrame:
    ''' 
    vol(): creates volcano plot
    
    Parameters:
    df (dataframe | str): pandas dataframe (or file path) from st.compare()
    FC (str): fold change column name (x-axis)
    pval (str): p-value column name (y-axis)
    stys (str, optional): style column name
    size (str | bool, optional): size column name (Default: -log10('pval'); specify False for no size)
    size_dims (tuple, optional): (minimum,maximum) values in size column (Default: None)
    z_col (str, optional): find rows with 'z_var' in this column for log2(FC) z-score normalization (Default: None)
    z_var (str, optional): use rows with 'z_var' for log2(FC) z-score normalization (Default: None)
    label (str, optional): label column name
    stys_order (list, optional): style column values order
    mark_order (list, optional): markers order for style column values order
    FC_threshold (float, optional): fold change threshold (Default: 1; log2(1)=0)
    pval_threshold (float, optional): p-value threshold (Default: 1; -log10(1)=0)
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
    legend_columnspacing (int, optional): space between columns (Default: -4; only for html plots)
    legend_handletextpad (float, optional): space between marker and text (Default: 0.5; only for html plots)
    legend_labelspacing (float, optional): vertical space between entries (Default: 0.5; only for html plots)
    legend_borderpad (float, optional): padding inside legend box (Default: 0.5; only for html plots)
    legend_handlelength (float, optional): marker length (Default: 0.5; only for html plots)
    legend_size_html_multiplier (float, optional): legend size multiplier for HTML plots (Default: 0.75)
    display_legend (bool, optional): display legend on plot (Default: True)
    display_labels (str | list, optional): display labels for values if label column specified (Options: 'FC & p-value', 'FC', 'p-value', 'NS', 'all', or [])
    display_lines (bool, optional): display lines for threshold (Default: False)
    display_axis (bool, optional): display x- and y-axis lines (Default: True)
    return_df (bool, optional): return dataframe (Default: True)
    dpi (int, optional): figure dpi (Default: 600 for non-HTML, 150 for HTML)
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    PDB_pt (str, optional): Path to PDB file for mol* visualization; only intended for fastq plots (Default: None)
    
    Dependencies: os, matplotlib, seaborn, & pandas
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)
    
    # Determine if we are saving to HTML (for interactive behavior)
    if file is not None:
        is_html = file.endswith('.html')
    else:
        is_html = False
    
    # Match title fontsize for html plots
    if is_html==True:
        x_axis_size=title_size
        y_axis_size=title_size
        x_ticks_size=title_size
        y_ticks_size=title_size
        legend_title_size=title_size
        legend_size=title_size*legend_size_html_multiplier
    
    # Prevent overwriting edms.fastq.vol()
    if f'log2({FC})_raw' not in df.columns:

        # Log transform data
        df[f'log2({FC})'] = [np.log10(FC_val)/np.log10(2) for FC_val in df[FC]]
        df[f'-log10({pval})'] = [-np.log10(pval_val) for pval_val in df[pval]]
        
        # Z-score normalization if specified
        if z_col is not None and z_var is not None:
            z_df = df[df[z_col]==z_var]
            mu = z_df[f'log2({FC})'].mean()
            sigma = z_df[f'log2({FC})'].std()
            df.rename(columns={f'log2({FC})': f'log2({FC})_raw'}, inplace=True)
            df[f'log2({FC})'] = [(val - mu) / sigma for val in df[f'log2({FC})_raw']]

    # Organize data by significance
    signif = []
    for (log2FC,log10P) in zip(df[f'log2({FC})'],df[f'-log10({pval})']):
        if (np.abs(log2FC)>=np.log10(FC_threshold)/np.log10(2))&(log10P>=-np.log10(pval_threshold)): signif.append(f'FC & p-value')
        elif (np.abs(log2FC)<np.log10(FC_threshold)/np.log10(2))&(log10P>=-np.log10(pval_threshold)): signif.append('p-value')
        elif (np.abs(log2FC)>=np.log10(FC_threshold)/np.log10(2))&(log10P<-np.log10(pval_threshold)): signif.append('FC')
        else: signif.append('NS')
    df['Significance']=signif
    
    # Organize data by 'pval' or specified 'size' column, typically input abundance
    sizes=(1,100)
    if size in [False,'False','false']: # No size
        size = None 

    else:
        if size is None: size = f'-log10({pval})' # default to pval

        if size is not None and size in df.columns:
            # Filter by size dimensions
            if size_dims is not None: 
                df = df[(df[size]>=size_dims[0])&(df[size]<=size_dims[1])]

            # Shared size normalization across all scatter calls so marker areas are consistent
            size_norm = None
            _vmin, _vmax = None, None
            if display_legend:
                _vmin = df[size].min()
                _vmax = df[size].max()
                # Guard against degenerate case where all values are equal
                if _vmin == _vmax:
                    _vmax = _vmin + 1e-12
                size_norm = mcolors.Normalize(vmin=_vmin, vmax=_vmax)
    
    # Set dimensions
    if x_axis_dims==(0,0): x_axis_dims=(min(df[f'log2({FC})']),max(df[f'log2({FC})']))
    if y_axis_dims==(0,0): y_axis_dims=(0,max(df[f'-log10({pval})']))

    # Generate figure
    fig, ax = plt.subplots(figsize=figsize)
    
    if vertical: # orientation
        # with significance boundraries
        if display_lines:
            plt.vlines(x=-np.log10(FC_threshold)/np.log10(2), ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
            plt.vlines(x=np.log10(FC_threshold)/np.log10(2), ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
            plt.hlines(y=-np.log10(pval_threshold), xmin=x_axis_dims[0], xmax=x_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)

        # with data
        if display_legend==False: size=None
        sns.scatterplot(
            data=df[df['Significance']!='FC & p-value'],
            x=f'log2({FC})', y=f'-log10({pval})',
            color=color, alpha=alpha, edgecolor=edgecol, 
            style=stys, style_order=stys_order if stys_order else None, markers=mark_order if mark_order else None, 
            size=size if display_legend else None, sizes=sizes, size_norm=size_norm,
            legend=False,
            ax=ax, **kwargs)
        sns.scatterplot(
            data=df[(df['Significance']=='FC & p-value')&(df[f'log2({FC})']<0)],
            x=f'log2({FC})', y=f'-log10({pval})',
            hue=f'log2({FC})', edgecolor=edgecol, palette='Blues_r', 
            style=stys, style_order=stys_order if stys_order else None, markers=mark_order if mark_order else None,
            size=size if display_legend else None, sizes=sizes, size_norm=size_norm,
            legend=False,
            ax=ax, **kwargs)
        sns.scatterplot(
            data=df[(df['Significance']=='FC & p-value')&(df[f'log2({FC})']>0)],
            x=f'log2({FC})', y=f'-log10({pval})',
            hue=f'log2({FC})', edgecolor=edgecol, palette='Reds', 
            style=stys, style_order=stys_order if stys_order else None, markers=mark_order if mark_order else None,
            size=size if display_legend else None, sizes=sizes, size_norm=size_norm, 
            legend=False,
            ax=ax, **kwargs)
        
        # with x- & y-axis lines
        if display_axis == True:
            ax.plot([x_axis_dims[0], x_axis_dims[1]], [0,0], color='black', linestyle='-', linewidth=1)
            ax.plot([0,0], [y_axis_dims[0], y_axis_dims[1]], color='black', linestyle='-', linewidth=1)

        
        # with legend
        if display_legend == True:
            if size_norm is not None and size is not None: # Add consistent size legend with 5 representative values
                if stys is not None and mark_order is not None and stys_order is not None: # Add stys legend too
                    legend_vals = np.linspace(_vmin, _vmax, len(mark_order))
                    for lv,mark,sty in zip(legend_vals,mark_order,stys_order):
                        plt.scatter([], [], s=np.interp(lv, [_vmin, _vmax], sizes), color=color, label=f'{lv:.2g}; {sty}', marker=mark)
                    if is_html:
                        if legend_bbox_to_anchor == (1,1): legend_bbox_to_anchor = (-0.1,-0.15)
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
                        plt.scatter([], [], s=np.interp(lv, [_vmin, _vmax], sizes), color=color, label=f'{lv:.2g}')
                    plt.legend(title=legend_title if legend_title!='' else size, 
                                title_fontsize=legend_title_size, fontsize=legend_size,
                                bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol)
        
        # with labels
        if label is not None and display_labels is not None:
            if display_labels in ['FC & p-value', 'FC', 'p-value', 'NS']:
                df_signif = df[df['Significance'] == display_labels]
            elif display_labels == 'all':
                df_signif = df
            elif isinstance(display_labels, list):
                df_signif = df[df[label].isin(display_labels)]
            else:
                df_signif = df[df['Significance'] == 'FC & p-value']
                print(f'Warning: defaulting to \"FC & p-value\" for label display due to invalid option for display_labels = {display_labels}')

            if is_html:
                # For HTML, show labels interactively as tooltips instead of static text
                pts = ax.scatter(
                    x=df_signif[f'log2({FC})'],
                    y=df_signif[f'-log10({pval})'],
                    s=20,
                    alpha=0
                )
                labels_list = df_signif[label].fillna("").astype(str).tolist()
                tooltip = SafeHTMLTooltip(pts, labels_list)
                clicker = ClickTooltip(pts, labels_list)
                mpld3.plugins.connect(fig, tooltip, clicker)
            else:
                # For static images, keep labels as always-visible text
                for i, l in enumerate(df_signif[label]):
                    plt.text(
                        x=df_signif.iloc[i][f'log2({FC})'],
                        y=df_signif.iloc[i][f'-log10({pval})'],
                        s=l
                    )
        
        # Set x axis
        if x_axis=='': x_axis=f'log2({FC})'
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
        if y_axis=='': y_axis=f'-log10({pval})'
        plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font,labelpad=y_axis_pad)

        if y_ticks==[]: plt.yticks(rotation=y_ticks_rot,fontfamily=y_ticks_font, fontsize=y_ticks_size)
        else: plt.yticks(ticks=y_ticks,labels=y_ticks,rotation=y_ticks_rot,fontfamily=y_ticks_font, fontsize=y_ticks_size)

    else: # Horizontal orientation
        # with significance boundraries
        if display_lines:
            plt.hlines(y=-np.log10(FC_threshold)/np.log10(2), xmin=y_axis_dims[0], xmax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
            plt.hlines(y=np.log10(FC_threshold)/np.log10(2), xmin=y_axis_dims[0], xmax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
            plt.vlines(x=-np.log10(pval_threshold), ymin=x_axis_dims[0], ymax=x_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)

        # with data
        if display_legend==False: size=None
        sns.scatterplot(
            data=df[df['Significance']!='FC & p-value'],
            y=f'log2({FC})', x=f'-log10({pval})',
            color=color, alpha=alpha, edgecolor=edgecol, 
            style=stys, style_order=stys_order if stys_order else None, markers=mark_order if mark_order else None,
            size=size if display_legend else None, sizes=sizes, size_norm=size_norm,
            legend=False,
            ax=ax, **kwargs)
        sns.scatterplot(
            data=df[(df['Significance']=='FC & p-value')&(df[f'log2({FC})']<0)],
            y=f'log2({FC})', x=f'-log10({pval})',
            hue=f'log2({FC})',edgecolor=edgecol, palette='Blues_r', 
            style=stys, style_order=stys_order if stys_order else None, markers=mark_order if mark_order else None,
            size=size if display_legend else None, sizes=sizes, size_norm=size_norm, 
            legend=False,
            ax=ax, **kwargs)
        sns.scatterplot(
            data=df[(df['Significance']=='FC & p-value')&(df[f'log2({FC})']>0)],
            y=f'log2({FC})', x=f'-log10({pval})',
            hue=f'log2({FC})', edgecolor=edgecol, palette='Reds', 
            style=stys, style_order=stys_order if stys_order else None, markers=mark_order if mark_order else None,
            size=size if display_legend else None, sizes=sizes, size_norm=size_norm, 
            legend=False,
            ax=ax, **kwargs)
        
        # with x- & y-axis lines
        if display_axis == True:
            ax.plot([y_axis_dims[0], y_axis_dims[1]], [0,0], color='black', linestyle='-', linewidth=1)
            ax.plot([0,0], [x_axis_dims[0], x_axis_dims[1]], color='black', linestyle='-', linewidth=1)

        # with legend
        if display_legend == True:
            if size_norm is not None and size is not None: # Add consistent size legend with 5 representative values
                if stys is not None and mark_order is not None and stys_order is not None: # Add stys legend too
                    legend_vals = np.linspace(_vmin, _vmax, len(mark_order))
                    for lv,mark,sty in zip(legend_vals,mark_order,stys_order):
                        plt.scatter([], [], s=np.interp(lv, [_vmin, _vmax], sizes), color=color, label=f'{lv:.2g}; {sty}', marker=mark)
                    if is_html:
                        if legend_bbox_to_anchor == (1,1): legend_bbox_to_anchor = (-0.1,-0.15)
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
                        plt.scatter([], [], s=np.interp(lv, [_vmin, _vmax], sizes), color=color, label=f'{lv:.2g}')
                    plt.legend(title=legend_title if legend_title!='' else size, 
                                title_fontsize=legend_title_size, fontsize=legend_size,
                                bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol)
        
        # with labels
        if label is not None:
            if display_labels in ['FC & p-value', 'FC', 'p-value', 'NS']:
                df_signif = df[df['Significance'] == display_labels]
            elif display_labels == 'all':
                df_signif = df
            elif isinstance(display_labels, list):
                df_signif = df[df[label].isin(display_labels)]
            else:
                df_signif = df[df['Significance'] == 'FC & p-value']
                print(f'Warning: defaulting to \"FC & p-value\" for label display due to invalid option for display_labels = {display_labels}')

            if is_html:
                # For HTML, show labels interactively as tooltips instead of static text
                pts = ax.scatter(
                    y=df_signif[f'log2({FC})'],
                    x=df_signif[f'-log10({pval})'],
                    s=20,
                    alpha=0
                )
                labels_list = df_signif[label].fillna("").astype(str).tolist()
                tooltip = SafeHTMLTooltip(pts, labels_list)
                clicker = ClickTooltip(pts, labels_list)
                mpld3.plugins.connect(fig, tooltip, clicker)
            else:
                # For static images, keep labels as always-visible text
                for i, l in enumerate(df_signif[label]):
                    plt.text(
                        y=df_signif.iloc[i][f'log2({FC})'],
                        x=df_signif.iloc[i][f'-log10({pval})'],
                        s=l
                    )
        
        # Set x axis
        if y_axis=='': y_axis=f'-log10({pval})'
        plt.xlabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font,labelpad=y_axis_pad)
        if y_ticks==[]: 
            if y_ticks_rot == 0: plt.xticks(rotation=y_ticks_rot,ha='center',va='top',fontfamily=y_ticks_font, fontsize=y_ticks_size)
            elif y_ticks_rot == 90: plt.xticks(rotation=y_ticks_rot,ha='right',va='center',fontfamily=y_ticks_font, fontsize=y_ticks_size)
            else: plt.xticks(rotation=y_ticks_rot,ha='right',fontfamily=y_ticks_font, fontsize=y_ticks_size)
        else: 
            if y_ticks_rot == 0: plt.xticks(ticks=y_ticks,labels=y_ticks,rotation=y_ticks_rot, ha='center',va='top',fontfamily=y_ticks_font, fontsize=y_ticks_size)
            elif y_ticks_rot == 90: plt.xticks(ticks=y_ticks,labels=y_ticks,rotation=y_ticks_rot,ha='right',va='center',fontfamily=y_ticks_font, fontsize=y_ticks_size)
            else: plt.xticks(ticks=y_ticks,labels=y_ticks,rotation=y_ticks_rot,ha='right',fontfamily=y_ticks_font, fontsize=y_ticks_size)
        
        # Set y axis
        if x_axis=='': x_axis=f'log2({FC})'
        plt.ylabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight,fontfamily=x_axis_font,labelpad=x_axis_pad)

        if x_ticks==[]: plt.yticks(rotation=x_ticks_rot,fontfamily=x_ticks_font, fontsize=x_ticks_size)
        else: plt.yticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot,fontfamily=x_ticks_font, fontsize=x_ticks_size)

    # Set title
    if title=='' and file is not None: 
        if space_capitalize: title=re_un_cap(".".join(file.split(".")[:-1]))
        else: ".".join(file.split(".")[:-1])
    plt.title(title, fontsize=title_size, fontweight=title_weight, family=title_font)

    # Save & show fig; return dataframe
    save_fig(file=file, dir=dir, fig=fig, dpi=dpi, PDB_pt=PDB_pt, icon='volcano')
    if show:
        ext = file.split('.')[-1].lower() if file is not None else ''
        if ext not in ('html', 'json'):
            plt.show()
        else:
            mpld3.show(fig)
    if return_df:
        return df

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