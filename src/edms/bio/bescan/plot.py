from .classes import *
from ...gen import io

def pwes_heatmap(
    df_scaled: pd.DataFrame | str,
    mask_on: bool = True,
    blackout: list[str] = ['X'],
    axis: AxisLabelOpts = AXIS,
    style: PWESHeatmapStyleOpts = PWESHEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    ):
    '''
    pwes_heatmap: Generates a heatmap for the scaled PWES matrix with options for masking, blackout regions, and styling.

    Parameters:
    - df_scaled (pd.DataFrame | str): Scaled DataFrame containing PWES values.
    - mask_on (bool): Whether to mask the lower triangle of the heatmap.
    - blackout (list[str]): List of prefixes for rows/columns to blackout (e.g., 'X' for non-protein residues).
    - axis (AxisLabelOpts): Axis labeling options.
    - style (PWESHeatmapStyleOpts): Styling options for the heatmap.
    - output (OutputOpts): Output options for saving and displaying the heatmap.
    '''
    # Get dataframe if path is provided #
    if isinstance(df_scaled, str):
        df_scaled = io.get(df_scaled, index_col=0)

    # SETUP #
    mask = np.zeros_like(df_scaled)
    if mask_on: mask[np.tril_indices_from(mask, k=-1)] = True # MASK FOR LOWER TRIANGLE #

    # GENERATE HEATMAP #
    fig, ax = plt.subplots(figsize=output.figsize)
    sns.heatmap(
        data=df_scaled,
        ax=ax,
        mask=mask,
        **style.heatmap_kws,
        )

    for edge, spine in ax.spines.items():
        spine.set_visible(style.spine_visible)
        spine.set_color(style.spine_color)

    # SET TICKS AND LABELS #
    ax.set_xlabel(axis.xlabel, fontproperties=arial_font6)
    ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)
    if axis.xticks: ax.set_xticklabels(axis.xticks)
    if axis.yticks: ax.set_yticklabels(axis.yticks)
    # TICK LABELS #
    for label in ax.get_xticklabels() + ax.get_yticklabels(): label.set_fontproperties(arial_font6)
    plt.subplots_adjust(wspace=style.spine_wspace, hspace=style.spine_hspace)

    # DRAW RECTANGLES #
    # HORIZONTAL #
    for i, idx in enumerate(df_scaled.index):
        for start in blackout:
            if idx.startswith(start):
                rect = patches.Rectangle(
                    (0, i), len(df_scaled.columns), 1,
                    linewidth=0, edgecolor=None, facecolor='black', alpha=1)
                ax.add_patch(rect)
    # VERTICAL #
    for j, col in enumerate(df_scaled.columns):
        for start in blackout:
            if col.startswith(start):
                rect = patches.Rectangle(
                    (j, 0), 1, len(df_scaled),
                    linewidth=0, edgecolor=None, facecolor='black', alpha=1)
                ax.add_patch(rect)

    # RASTERIZE #
    if output.rasterize:
        for coll in ax.collections:
            coll.set_rasterized(True)

    plt.tight_layout()
    if output.save:
        out_path = Path(f"{str(output.path)}_heatmap.{output.out_type}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, format=output.out_type,
                    dpi=output.dpi, transparent=output.transparent
                    )
    # SHOW #
    if output.show: plt.show()
    plt.close(fig)

def pwes_clustermap(
    df_scaled, link, df_clusters,
    color_list=[],

    axis: AxisLabelOpts = AXIS,
    style: PWESHeatmapStyleOpts = PWESHEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    ):

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
        figsize=output.figsize,
        **style.clustered_heatmap_kws,
        )
    ax = fig.ax_heatmap

    # AXIS SETTINGS #
    fig.ax_heatmap.collections[0].set_rasterized(output.rasterize)
    fig.cax.set_visible(False)
    ax.set_aspect('equal')
    fig.ax_col_dendrogram.set_visible(False)

    # DRAW LINES BETWEEN CLUSTERS #
    temp = 0
    clusterlist = df_clusters['cl_new'].value_counts().sort_index()
    for i in clusterlist.iloc[:-1]:
        temp += i
        ax.axhline(y=temp, **style.line_kws)

    for edge, spine in ax.spines.items():
        spine.set_visible(style.spine_visible)
        spine.set_color(style.spine_color)

    plt.tight_layout()

    # SAVE #
    if output.save:
        out_path = Path(f"{str(output.path)}_cluster_heatmap.{output.out_type}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, format=output.out_type,
                    dpi=output.dpi, transparent=output.transparent
                    )
    # SHOW #
    if output.show: plt.show()
    plt.close()

def pwes_heatmap_colorbar(
    orientation='vertical', # 'v' or 'h'

    axis: AxisLabelOpts = AXIS,
    style: PWESHeatmapStyleOpts = PWESHEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    fig, ax = plt.subplots(figsize=output.figsize)

    norm = mpl.colors.Normalize(**style.colorbar_normalize_kws)
    cbar = mpl.colorbar.ColorbarBase(
        ax, norm=norm, orientation=orientation,
        **style.colorbar_base_kws)

    for label in (ax.get_yticklabels() if orientation == 'vertical' else ax.get_xticklabels()):
        label.set_fontproperties(arial_font6)

    ax.tick_params(**axis.tick_kws)
    cbar.outline.set_linewidth(style.colorbar_linewidth)

    if output.save:
        out_path = Path(f"{str(output.path)}.{output.out_type}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, format=output.out_type,
                    dpi=output.dpi, transparent=output.transparent
                    )
    # SHOW #
    if output.show: plt.show()
    plt.close(fig)

# PWES Cluster Mapping #

def tuple_to_hex(rgb_tuple):
    if not (isinstance(rgb_tuple, tuple) and len(rgb_tuple) == 3 and
            all(isinstance(val, float) and 0 <= val <= 1 for val in rgb_tuple)):
        return None

    r, g, b = rgb_tuple
    r, g, b = int(r * 255), int(g * 255), int(b * 255)
    return f"{r:02x}{g:02x}{b:02x}".upper()

def generate_pymol_script(
    output_filename,
    pdb_filename,
    residue_dict,
    colors, 
    ):

    """
    Generates a PyMOL script to display selected residues as separate objects.
    """

    with open(output_filename, 'w') as pymol_script:
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
