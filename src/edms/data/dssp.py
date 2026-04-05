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
- ssa_str(): Convert DSSP secondary structure codes for a chain into a string format compatible with secstructartist (GitHub).
- dssp_pdb_id(): Extract the PDB ID from a DSSP file header, if present.

[Plot functions]
- plot_ss_color_key(): Plot a color key for DSSP secondary structure assignments.
- plot_ssa(): Plot the secondary structure assignment for a given DSSP file and chain.

[Pymol functions]
- pymol_color_defs_from_ss_map(): Convert ss_map matplotlib colors into PyMOL set_color commands.
- pymol_ssa(): Write a PyMOL script that colors secondary-structure elements from a DSSP file.
'''
# Import packages
import re
import pandas as pd
import subprocess
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgb
from secstructartist.artists import (
        SecStructArtist, 
        ElementArtist, 
        primitives as prim,
    )
from pathlib import Path
from ..utils import mkdir
from ..gen import plot as p

# DSSP code descriptions
ss_map = {
    # Helices
    "H": ["α-helix", "turquoise"], # i+4 helix
    "G": ["3_10-helix", "mediumblue"], # i+3 helix
    "I": ["π-helix", "slateblue"], # i+5 helix
    "P": ["κ-helix", "forestgreen"], # polyproline helix
    # Strands
    "B": ["β-bridge", "orange"], # single residue β-strand
    "E": ["β-strand", "gold"], # extended β-strand
    # Flexible
    "T": ["turn", "salmon"], # hydrogen-bonded turn
    "S": ["bend", "sienna"], # bend (non-hydrogen-bonded turn)
    " ": ["coil", "darkgray"], # loop or irregular structure
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

def ssa_str(dssp_file: str, chain_id: str,) -> str:
    """
    ssa_str(): Convert DSSP secondary structure codes for a chain into a string format compatible with secstructartist (GitHub).

    Parameters:
    dssp_file (str): Path to legacy DSSP output file.
    chain_id (str): Chain identifier to extract, e.g. "A".
    """
    lines = Path(dssp_file).read_text().splitlines()

    # DSSP data starts after the header line beginning with "  #"
    start = None
    for i, line in enumerate(lines):
        if line.startswith("  #"):
            start = i + 1
            break
    if start is None:
        raise ValueError("Could not find DSSP residue table header line starting with '  #'.")

    out = []

    for line in lines[start:]:
        if not line.strip():
            continue
        if len(line) < 17:
            continue

        # In legacy DSSP format:
        # chain id is around column 12 (0-based index 11)
        # SS code is around column 17 (0-based index 16)
        # These positions are standard for classic DSSP text output.
        this_chain = line[11]
        ss_code = line[16]

        # Skip rows from other chains
        if this_chain != chain_id:
            continue

        out.append(ss_code)

    return "".join(out)

def dssp_pdb_id(dssp_file: str) -> str | None:
    """
    dssp_pdb_id(): Extract the PDB ID from a DSSP file header, if present.

    Parameters:
    dssp_file (str): Path to the DSSP file.
    """
    with open(dssp_file, "r") as f:
        for line in f:
            if line.startswith("HEADER"):
                # Try fixed-width PDB-style position first.
                # In classic PDB HEADER records, the ID code is near the end.
                fixed = line[62:66].strip().upper() if len(line) >= 66 else ""
                if re.fullmatch(r"[A-Z0-9]{4}", fixed):
                    return fixed

                # Fallback: search the HEADER line for 4-char alphanumeric tokens.
                tokens = re.findall(r"\b[A-Za-z0-9]{4}\b", line.upper())
                if tokens:
                    return tokens[-1]

                return None

    return None

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
    
    plot_ss_color_key

def plot_ssa(dssp_file: str, chain_id: str, artist: SecStructArtist | str = None, figsize: tuple = (8, 0.5),
            dir: str = None, file: str = None, dpi: int = 1200, transparent: bool = True, show: bool = True, **kwargs):
    """
    plot_ssa(): Plot the secondary structure assignment for a given DSSP file and chain.

    Parameters:
    dssp_file (str): Path to the DSSP file.
    chain_id (str): Chain identifier to extract, e.g. "A".
    artist (SecStructArtist | str): Optional secstructartist ElementArtist to use for plotting; if None, edms custom style will be used.
    dir (str, optional): Directory to save the figure; if None, the figure is not saved.
    file (str, optional): Filename to save the figure; if None, the figure is not saved.
    dpi (int, optional): Resolution for saving the figure.
    transparent (bool, optional): Whether to save the figure with a transparent background.
    show (bool, optional): Whether to display the plot (default: True).
    """
    # Custom artist (default if not provided)
    ARTIST_NAME = 'dssp-custom'

    ### Definition of the Element artists
    H = ElementArtist([
            prim.HelixPrimitive(
                fillcolor='turquoise',
                linecolor='#111111',
                height_scalar=.9, 
                ribbon_period = 3.6,
                fill_inner_ribbon = False,
            ),
            prim.HelixPrimitive(
                fillcolor="paleturquoise",
                linecolor='#111111',
                height_scalar=.9, 
                ribbon_period = 3.6,
                fill_inner_ribbon = True,
                zorder_offset = -.5
            )
        ], r'$\alpha$-Helix'
    )
    G = ElementArtist([
            prim.HelixPrimitive(
                fillcolor='mediumblue',
                linecolor='#111111',
                height_scalar=.7, 
                ribbon_period = 3.0,
                fill_inner_ribbon = False,
            ),
            prim.HelixPrimitive(
                fillcolor="royalblue",
                linecolor='#111111',
                height_scalar=.7, 
                ribbon_period = 3.0,
                fill_inner_ribbon = True,
                zorder_offset = -.5
            )
        ], r'$3_{10}$-Helix'
    )
    P = ElementArtist([
            prim.HelixPrimitive(
                fillcolor='forestgreen',
                linecolor='#111111',
                height_scalar=.5, 
                ribbon_period = 3.0,
                fill_inner_ribbon = False,
            ),
            prim.HelixPrimitive(
                fillcolor="palegreen",
                linecolor='#111111',
                height_scalar=.5, 
                ribbon_period = 3.0,
                fill_inner_ribbon = True,
                zorder_offset = -.5
            )
        ], r'$\kappa$-Helix'
    )
    I = ElementArtist([
            prim.HelixPrimitive(
                fillcolor='slateblue',
                linecolor='#111111',
                height_scalar=1., 
                ribbon_period = 4.1,
                fill_inner_ribbon = False,
            ),
                prim.HelixPrimitive(
                fillcolor="mediumslateblue",
                linecolor='#111111',
                height_scalar=1., 
                ribbon_period = 4.1,
                fill_inner_ribbon = True,
                zorder_offset = -.5
            )
        ], r'$\pi$-Helix'
    )
    B = ElementArtist([
            prim.ArrowPrimitive(
                fillcolor='orange',
                linecolor='#111111', 
                height_scalar = .65,
                height_scalar2 = .4,
                arrow_tip_length =3
            )
        ], r'$\beta$-Bridge'
    )
    E = ElementArtist([
            prim.ArrowPrimitive(
                fillcolor='gold',
                linecolor='#111111', 
                height_scalar = .65,
                height_scalar2 = .4,
                arrow_tip_length = 3
            )
        ], r'Extended $\beta$-Strand'
    )
    S = ElementArtist([
            prim.LinePrimitive(
                linecolor = "sienna", 
                linewidth_scalar = 5
            )
        ], 'Bend'
    )
    T = ElementArtist([
            prim.LinePrimitive(
                linecolor="salmon", 
                linewidth_scalar = 5
            )
        ], 'H-bonded turn'
    )
    C = ElementArtist([
            prim.LinePrimitive(
                linecolor="black", 
                linewidth_scalar = 1
            )
        ], 'Loop'
    )

    ### Artist setup
    dssp_mapping = {
        'H': H, 'h': H, 
        'G': G, 'g': G, 
        'P': P, 'p': P, 
        'I': I, 'i': I,
        'E': E, 'e': E, 
        'B': B, 'b': B, 
        'T': T, 't': T, 
        'S': S, 's': S, 
        'C': C, 'c': C, ' ': C, '-': C,
    }

    # Use custom artist if not provided
    if artist is None:
        artist = SecStructArtist(dssp_mapping, linewidth=.5, zorder=5)
    
    # Plotting
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': 'secstruct'})
    ax.draw_secondary_structure(
        secstruct=ssa_str(dssp_file=dssp_file, chain_id=chain_id),
        artist=artist
    )

    # Remove the axis
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_axis_off()

    # Legend (only viewable with jupyter notebook backend)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.6), ncol=3)
    
    plt.tight_layout()
    if dir is not None and file is not None:
        p.save_fig(dir=dir, file=file, dpi=dpi, transparent=transparent)
    if show:
        plt.show()

# Pymol functions
def pymol_color_defs_from_ss_map(ss_map: dict) -> list[str]:
    """
    pymol_color_defs_from_ss_map(): Convert ss_map matplotlib colors into PyMOL set_color commands.

    Parameters:
    ss_map (dict): Mapping of DSSP codes to descriptions and matplotlib color names.
    """
    lines = []
    for code, (_, mpl_color) in ss_map.items():
        pymol_name = f"dssp_{code}" if code != " " else "dssp_coil"
        r, g, b = to_rgb(mpl_color)
        lines.append(f"set_color {pymol_name}, [{r:.4f}, {g:.4f}, {b:.4f}]")
    return lines

def pymol_ssa(
    dssp_file: str,
    chain_id: str,
    dir: str,
    file: str,
    pdb_id_or_filename: str = None,
    object_name: str = "prot",
    base_color: str = "white",
    unknown_color: str = "white",
    execute: bool = False,
):
    """
    pymol_ssa(): Write a PyMOL script that colors secondary-structure elements from a DSSP file.

    Parameters:
    dssp_file (str): Path to the DSSP file.
    chain_id (str): Chain identifier to extract, e.g. "A".
    dir (str): Directory to save the PyMOL script.
    file (str): Output PyMOL script filename, e.g. 'ssa.pml'.
    pdb_id_or_filename (str): PDB ID or path to the structure file to load in PyMOL.
    object_name (str, optional): Name of the PyMOL object.
    base_color (str, optional): Base color applied before DSSP coloring.
    unknown_color (str, optional): Color for unknown DSSP codes.
    execute (bool, optional): Whether to execute the generated PyMOL script immediately after writing it.
    """
    mkdir(dir)
    out_file = os.path.join(dir, file)

    df = parse_segments(
        dssp_file=dssp_file,
        chain_id=chain_id,
        unknown_color=unknown_color,
    ).copy()

    if df.empty:
        raise ValueError(f"No DSSP segments found for chain '{chain_id}'.")

    def _sanitize_name(text: str) -> str:
        text = str(text)
        return (
            text.replace("α", "alpha")
                .replace("β", "beta")
                .replace("π", "pi")
                .replace("κ", "kappa")
                .replace(" ", "_")
                .replace("-", "_")
                .replace("/", "_")
                .replace("(", "")
                .replace(")", "")
                .replace(".", "")
                .replace(",", "")
                .replace("'", "")
                .replace("#","")
                .lower()
        )

    def _pymol_color_name(ss_code: str) -> str:
        return f"dssp_{ss_code}" if ss_code != " " else "dssp_coil"

    flexible_codes = {"T", "S", " "}

    with open(out_file, "w") as pymol_script:
        if pdb_id_or_filename is None:
            pdb_id_or_filename = dssp_pdb_id(dssp_file = dssp_file)
            if pdb_id_or_filename is None:
                print("Warning: Could not determine PDB ID from DSSP file. Please provide a structure file path or valid PDB ID to load in PyMOL.")
            pymol_script.write(f"fetch {pdb_id_or_filename}, {object_name}\n")
        elif len(pdb_id_or_filename) == 4 and pdb_id_or_filename.isalnum():
            pymol_script.write(f"fetch {pdb_id_or_filename}, {object_name}\n")
        else:
            pymol_script.write(f'load "{os.path.abspath(pdb_id_or_filename)}", {object_name}\n')
        pymol_script.write(f"hide everything, {object_name}\n")
        pymol_script.write(f"show cartoon, {object_name}\n")
        pymol_script.write(f"color {base_color}, {object_name}\n")
        pymol_script.write("bg_color white\n")
        pymol_script.write("set cartoon_fancy_helices, 1\n")
        pymol_script.write("set cartoon_flat_sheets, 1\n")
        pymol_script.write("set cartoon_side_chain_helper, 0\n\n")

        # Define PyMOL colors from matplotlib color names in ss_map
        pymol_script.write("# DSSP colors converted from matplotlib names\n")
        for line in pymol_color_defs_from_ss_map(ss_map):
            pymol_script.write(line + "\n")
        pymol_script.write("\n")

        # Optional fallback color for unknown DSSP codes
        r, g, b = to_rgb(unknown_color)
        pymol_script.write(f"set_color dssp_unknown, [{r:.4f}, {g:.4f}, {b:.4f}]\n\n")

        # Color DSSP segments directly on the original object
        for row in df.itertuples(index=False):
            start = int(row.start)
            end = int(row.end)
            ss_code = row.ss_code

            selection = f"{object_name} and chain {chain_id} and resi {start}-{end}"
            pymol_color_name = _pymol_color_name(ss_code) if ss_code in ss_map else "dssp_unknown"

            pymol_script.write(f"color {pymol_color_name}, ({selection})\n")

            pymol_script.write("\n")

        # Optional selections by DSSP code
        pymol_script.write("# Group residues by DSSP code\n")
        for code in df["ss_code"].drop_duplicates():
            code_df = df[df["ss_code"] == code]
            if code_df.empty:
                continue

            sel_parts = [
                f"(chain {chain_id} and resi {int(r.start)}-{int(r.end)})"
                for r in code_df.itertuples(index=False)
            ]

            code_name = "coil" if code == " " else _sanitize_name(ss_map.get(code, [code])[0])
            pymol_script.write(f"select dssp_{code_name}, " + " or ".join(sel_parts) + "\n")

        pymol_script.write("\n")
        pymol_script.write(f"orient {object_name} and chain {chain_id}\n")
        pymol_script.write(f"zoom {object_name} and chain {chain_id}\n")
        pymol_script.write("deselect\n")
    
    if execute:
        subprocess.run(f"pymol {out_file}", shell=True)