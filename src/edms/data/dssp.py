'''
Module: dssp.py
Author: Marc Zepeda
Created: 2026-03-29
Description: Define Secondary Structure of Proteins

[DSSP code descriptions]
- ss_map: Mapping of DSSP secondary structure codes to descriptions and colors.

[DSSP functions]
- retrieve(): Retrieve DSSP data for a given PDB ID.
- parse_segments(): Read a DSSP file, isolate one chain, and return continuous DSSP segments.

[Plot functions]
- plot_ss_color_key(): Plot a color key for DSSP secondary structure assignments.
'''
# Import packages
import pandas as pd
import subprocess
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from ..utils import mkdir
from ..gen import plot as p

# DSSP code descriptions
ss_map = {
    "H": ["α-helix", "turquoise"], # Helix
    "G": ["3_10-helix", "mediumblue"],
    "I": ["π-helix", "slateblue"],
    "P": ["polyproline-helix", "forestgreen"],
    "B": ["β-bridge", "orange"], # Strand
    "E": ["β-strand", "gold"],
    "T": ["turn", "salmon"], # Flexible
    "S": ["bend", "sienna"],
    " ": ["coil", "darkgray"],
}

# DSSP functions
def retrieve(id_or_file: str, dir: str=None) -> pd.DataFrame:
    """
    retrieve(): Retrieve DSSP data for a given PDB ID.

    Parameters:
    id_or_file (str): 4-character PDB accession (e.g. "1CRN", case-insensitive).
    dir (str, optional): Save file to specified directory; otherwise, save to configuration directory (~/.config/edms/DSSP/).
    
    """
    # Make directories if needed
    if dir is None:
        dir = os.path.expanduser('~/.config/edms/DSSP')
    mkdir(os.path.expanduser('~/.config/edms/DSSP/'))
    
    # Load PDB structure data
    if os.path.isfile(id_or_file): # from filepath
        command = f"mkdssp {id_or_file} {dir}/{os.path.basename(id_or_file).split('.')[0]}.dssp"
        print(f"Running command: {command}")
        subprocess.run(f"{command}", shell=True, cwd='.')
    else:
        try: # from config
            for id_or_file_config in os.listdir(os.path.expanduser('~/.config/edms/PDB/')):
                if id_or_file.lower() in id_or_file_config.lower() and (id_or_file_config.endswith('.cif') or id_or_file_config.endswith('.mmCIF')):
                    command = f"mkdssp {os.path.expanduser('~/.config/edms/PDB')}/{id_or_file_config} {dir}/{id_or_file}.dssp"
                    print(f"Running command: {command}")
                    subprocess.run(f"{command}", shell=True, cwd='.')
                    break

        except:
            raise FileNotFoundError(f"PDB file with .cif or .mmCIF extension was not found: {id_or_file}.\nPlease provide a valid filename or PDB id (if saved to {os.path.expanduser('~/.config/edms/PDB/')}) or file path for PDB structure file")


def parse_segments(dssp_file: str, chain_id: str, unknown_color: str='white') -> pd.DataFrame:
    """
    parse_segments(): Read a DSSP file, isolate one chain, and return continuous DSSP segments.

    Parameters:
    dssp_file (str): Path to the DSSP file.
    chain_id (str): Chain identifier to isolate.
    unknown_color (str): Color to use for unknown secondary structure codes (default: 'white')
    """
    records = []
    in_data = False

    with open(dssp_file, "r") as f:
        for line in f:
            if line.startswith("  #  RESIDUE AA STRUCTURE"):
                in_data = True
                continue

            if not in_data or not line.strip():
                continue

            try:
                resnum_str = line[5:10].strip()
                chain = line[11].strip()
                aa = line[13].strip()
                ss = line[16].strip()
            except IndexError:
                continue

            if not resnum_str:
                continue

            if aa == "!":
                continue

            if ss == "":
                ss = " "

            try:
                resnum = int(resnum_str)
            except ValueError:
                continue

            records.append(
                {
                    "resnum": resnum,
                    "chain": chain,
                    "aa": aa,
                    "ss_code": ss,
                }
            )

    df = pd.DataFrame(records)

    if df.empty:
        raise ValueError("No DSSP residue records were parsed.")

    df = df[df["chain"] == chain_id].copy()

    if df.empty:
        raise ValueError(f"No residues found for chain '{chain_id}'.")

    df = df.sort_values("resnum").reset_index(drop=True)

    df["new_segment"] = (
        (df["ss_code"] != df["ss_code"].shift(1))
        | (df["resnum"] != df["resnum"].shift(1) + 1)
    )

    df["segment_id"] = df["new_segment"].cumsum()

    out = (
        df.groupby(["segment_id", "ss_code"], as_index=False)
        .agg(start=("resnum", "min"), end=("resnum", "max"))
        .copy()
    )

    out["chain"] = chain_id

    # add color column
    out["ss_color"] = out["ss_code"].map(ss_map).apply(lambda x: x[1] if isinstance(x, list) else unknown_color)

    # add description column
    out["ss_description"] = out["ss_code"].map(ss_map).apply(lambda x: x[0] if isinstance(x, list) else x)

    # count per structure type
    out["ss_count"] = out.groupby("ss_code").cumcount() + 1

    # add name column
    out["ss_name"] = out["ss_description"] + " #" + out["ss_count"].astype(str)

    out = out[
        ["chain", "ss_code", "ss_color", "ss_description", "ss_count", 'ss_name', "start", "end"]
    ]

    return out

# Plot functions
def plot_ss_color_key(
    ncols: int = 3,
    figsize: tuple = (8, 3),
    title: str = "2° Structure",
    title_size: int = 12,
    title_weight: str = 'bold',
    fontsize: int = 12,
    swatch_width: float = 0.35,
    swatch_height: float = 0.35,
    alpha: float = 0.15,
    dir: str = None,
    file: str = None,
    dpi: int = 1200,
    transparent: bool = True,
    show: bool = True
):
    """
    plot_ss_color_key(): Plot a color key for DSSP secondary structure assignments.

    Parameters:
    ncols (int, optional): Number of columns in the key.
    figsize (tuple, optional): Figure size.
    title (str, optional): Figure title.
    title_size (int, optional): Title font size.
    fontsize (int, optional): Label font size.
    swatch_width (float, optional): Width of color box.
    swatch_height (float, optional): Height of color box.
    alpha (float, optional): Transparency of color boxes.
    dir (str, optional): Directory to save the figure; if None, the figure is not saved.
    file (str, optional): Filename to save the figure; if None, the figure is not saved.
    dpi (int, optional): Resolution for saving the figure.
    transparent (bool, optional): Whether to save the figure with a transparent background.
    show (bool, optional): Whether to display the figure.
    """
    ss_colors = {}
    ss_labels = {}
    for code, (desc, color) in ss_map.items():
        ss_colors[code] = color
        ss_labels[code] = desc

    codes = list(ss_colors.keys())
    n_items = len(codes)
    nrows = (n_items + ncols - 1) // ncols

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, ncols * 2.8)
    ax.set_ylim(0, nrows * 1.0)
    ax.axis("off")

    for i, code in enumerate(codes):
        row = i // ncols
        col = i % ncols

        x = col * 2.8
        y = nrows - row - 0.7

        rect = Rectangle(
            (x, y),
            swatch_width,
            swatch_height,
            facecolor=ss_colors[code],
            edgecolor="black",
            linewidth=0.8,
            alpha=alpha
        )
        ax.add_patch(rect)

        label = f"{ss_labels.get(code, code)}"
        ax.text(
            x + swatch_width + 0.15,
            y + swatch_height / 2,
            label,
            va="center",
            ha="left",
            fontsize=fontsize,
        )

    ax.set_title(title, fontsize=title_size, fontweight=title_weight)
    
    plt.tight_layout()
    if dir is not None and file is not None:
        p.save_fig(dir=dir, file=file, dpi=dpi, transparent=transparent)
    if show:
        plt.show()