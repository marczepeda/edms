''' 
Module: genbank.py
Author: Marc Zepeda
Created: 2024-12-15
Description: Genbank visualization

Usage:
[Visualize]
- viewer(): visualize genbank file
'''
# Import packages
from dna_features_viewer import GraphicFeature, GraphicRecord
from Bio import SeqIO
import matplotlib.pyplot as plt
import os
from ..gen import io
from ..gen import plot as p

# Visualize
def viewer(pt: str, feature_colors: dict=None, exclude: list=[], region:tuple=None,
           file: str=None, dir: str=None, fig_width: int=10,
           title: str='',title_size: int=18, title_weight: str='bold',
           x_axis: str='bp', x_axis_size: int=12, x_axis_weight: str='bold', xticks: list=[],
           legend_title: str=None, legend_title_size: int=12,legend_bbox_to_anchor: tuple=(0.5,-0.25),
           legend_loc: str='upper center',legend_frame_on: bool=True,legend_ncol: int=4,
           return_features: bool=False, show: bool=True) -> list[tuple[str, GraphicFeature]]:
    '''
    viewer(): visualize genbank file

    Parameters:
    pt (str): Genbank file path
    feature_colors (dict, optional): dictionary of features and corresponding hexadecimal colors (Default: None; Ex: {'feature':'#000000'})
    exclude (list, optional): list of excluded features/annotations (Default: [])
    region (tuple, optional): (start, end) coordinates for highlighed region (Default: None)
    file (str, optional): save plot to filename (Default: None)
    dir (str, optional): save plot to directory (Default: None)
    fig_width (int | float, optional): figure width (Default: 10)
    title (str, optional): plot title (Default: '')
    title_size (int, optional): plot title font size (Default: 18)
    title_weight (str, optional): plot title bold, italics, etc. (Default: 'bold')
    x_axis (str, optional): x-axis name (Default: 'bp')
    x_axis_size (int, optional): x-axis name font size (Default: 12)
    x_axis_weight (str, optional): x-axis name bold, italics, etc. (Default: 'bold')
    xticks (list, optional): x-axis tick values (Default: [])
    legend_title (str, optional): legend title (Default: None)
    legend_title_size (str, optional): legend title font size (Default: 12)
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor (Default: (0.5,-0.25))
    legend_loc (str, optional): legend location (Default: 'upper center')
    legend_frame_on (bool, optional): display frame around legend (Default: True)
    legend_ncol (tuple, optional): # of columns (Default: 4)
    return_features (bool, optional): return features list? (Default: False)
    show (bool, optional): show plot (Default: True)

    Dependencies: Bio,dna_features_viewer,matplotlib
    '''
    # Parse GenBank file
    record = SeqIO.read(pt, "genbank")
    
    # Define a color scheme for feature types
    color_scheme = {
        "3' UTR": "#F5F5F5",            # White
        "5' UTR": "#F5F5F5",            # White
        "CDS": "#FFFF00",               # Yellow
        "DNA Probe": "#FF0000",         # Red
        "enhancer": "#FFC800",          # Yellow
        "enzyme": "#0000FF",            # Blue
        "exon": "#646464",              # Black
        "gene": "#00B200",              # Green
        "intron": "#F5F5F5",            # White
        "LTR": "#E7960A",               # Orange
        "mat_peptide": "#005A28",       # Dark green
        "misc_feature": "#AAAAAA",      # Gray
        "misc_difference": "#2E96FF",   # Light blue
        "misc_RNA": "#FF0064",          # Light red
        "modified_base": "FF6400",      # Orange
        "motif": "#96008C",             # Maroon
        "mRNA": "#BE3232",              # Brown-red
        "mutation": "#00009D",          # Dark blue
        "ncRNA": "#BE6464",             # Clay
        "orf": "#FF9900",               # Orange
        "rep_origin": "#50B3FF",        # Light blue
        "overhang": "#95E2FF",          # Lighter blue
        "polyA_signal": "#00CCCC",      # Turquoise
        "polyA_site": "#00FFFF",        # Cyan
        "polylinker": "#7100A8",        # Purple
        "polymorphism": "#FFB400",      # Orange
        "primer_bind": "#007C00",       # Dark green
        "promoter": "#BAFF00",          # Neon green
        "regulatory": "#0EA089",        # Dark turquoise
        "repeat_region": "#E7960A",     # Orange
        "rRNA": "#FF0000",              # Red
        "sig_peptide": "#FF73FF",       # Pink
        "source": "#0000FF",            # Deep Ocean Blue
        "terminator": "#FF5300",        # Orange
        "tmRNA": "#BE3232",             # Brown
        "Unsure": "#FFB400",            # Yellow-orange
        "Zinc Finger": "#DCAA8C"        # Tan
    }
    if feature_colors is not None: 
        for feature,color in feature_colors.items():
            color_scheme[feature]=color

    # Extract features (i.e., annotations) and assign colors based on their type
    features = []
    for feature in record.features: # Iterate through features
        if region is None: # Full sequence
            if feature.type not in exclude: # Check if features should be excluded
                features.append((feature.type,
                                GraphicFeature(start=int(feature.location.start), # Get feature information
                                               end=int(feature.location.end), 
                                               strand=feature.location.strand, 
                                               label=feature.qualifiers['label'][0], 
                                               color=color_scheme.get(feature.type, "#cccccc")))) # Default color (gray)
        else: # Specified region
            if feature.type not in exclude: # Check if features should be excluded
                if ((int(feature.location.start)>=region[0])&(int(feature.location.start)<=region[1]))|((int(feature.location.end)>=region[0])&(int(feature.location.end)<=region[1])): # Check if start and/or end of the feature is within the region
                    features.append((feature.type,
                                    GraphicFeature(start=int(feature.location.start), # Get feature information
                                                   end=int(feature.location.end), 
                                                   strand=feature.location.strand, 
                                                   label=feature.qualifiers['label'][0], 
                                                   color=color_scheme.get(feature.type, "#cccccc")))) # Default color (gray)

    # Create a GraphicRecord for visualization
    graphic_record = GraphicRecord(sequence_length=len(record.seq), features=[feature[1] for feature in features])
    graphic_record.plot(figure_width=fig_width)

    # Set title and x-axis
    if title=='' and file is not None: title=p.re_un_cap(".".join(file.split(".")[:-1]))
    elif title=='': title=p.re_un_cap(".".join(pt.split("/")[-1].split(".")[:-1]))
    plt.title(title, fontsize=title_size, fontweight=title_weight)
    plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight)
    plt.tick_params(axis='x', colors='black')
    if xticks!=[]: plt.xticks(ticks=xticks,labels=xticks)
    if region is not None: # Highlight specified region
        plt.xlim(region[0],region[1])

    # Add a legend
    if legend_title is not None:
        feature_types = sorted({(feature[0], color_scheme.get(feature[0], "#cccccc")) for feature in features}) # Get features and associated colors
        legend_elements = [plt.Line2D([0], [0], color=color, lw=4, label=label) for label,color in feature_types] # Make legend
        plt.legend(handles=legend_elements, bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, fontsize=legend_title_size, frameon=legend_frame_on, title=legend_title, ncol=legend_ncol)


    # Save, show fig, & return features
    if file is not None and dir is not None:
        io.mkdir(dir) # Make output directory if it does not exist
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight', format=f'{file.split(".")[-1]}')
    if show: plt.show()
    if return_features: return features