'''
src/edms/data/cli.py             Command Line Interface for EDMS Database module
├── __init__.py                 Initializer
├── cosmic.py                   COSMIC module
├── cvar.py                     ClinVar module
├── cBioPortal.py               cBioPortal module
├── uniprot.py                  Uniprot module
├── pdb.py                      PDB module
└── dssp.py                     DSSP module
'''
import argparse
import datetime
import sys # might use later
from rich import print as rprint # might use later

from . import cosmic as co, cvar, cBioPortal as cBP, uniprot, pdb, dssp
from ..utils import parse_tuple_int, parse_tuple_float # might use later

def add_subparser(subparsers, formatter_class=None):
    """
    add_subparser(): Attach all dat-related subparsers to the top-level CLI.

    Parameters:
    subparsers (argparse._SubParsersAction): The subparsers object to attach the dat subparsers to.
    formatter_class (type, optional): The formatter class to use for the subparsers.
    """
    if formatter_class is None:
        # fall back to basic formatter to avoid circular imports
        formatter_class = argparse.HelpFormatter

    '''
    edms.data.cosmic:
    - mutations(): returns COSMIC mutations dataframe for a given gene
    - cds_group(): plot COSMIC mutations histogram with CDS regions highlighted in different colors
    - priority_muts(): returns the shared sequences library dataframe with priority mutations
    - priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    - editor_mutations(): returns and plots editor accessible COSMIC mutations
    '''
    parser_cosmic = subparsers.add_parser("cosmic", help="COSMIC Database", description="COSMIC Database", formatter_class=formatter_class)
    subparsers_cosmic = parser_cosmic.add_subparsers()

    # mutations(): returns COSMIC mutations dataframe for a given gene
    parser_cosmic_mutations = subparsers_cosmic.add_parser("mutations", help="Extract COSMIC mutations", description="Extract COSMIC mutations", formatter_class=formatter_class)

    parser_cosmic_mutations.add_argument("-d", "--df", type=str, help="Input file path", required=True)

    parser_cosmic_mutations.add_argument("-o", "--dir", type=str, help="Output directory", default='../out')
    parser_cosmic_mutations.add_argument("-f", "--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_cosmic_mutations.csv')

    parser_cosmic_mutations.set_defaults(func=co.mutations)
    
    # cds_group(): plot COSMIC mutations histogram with CDS regions highlighted in different colors
    parser_cds_group = subparsers_cosmic.add_parser("cds_group", help="Plot COSMIC mutation histogram with CDS regions highlighted", description="Plot COSMIC mutation histogram with CDS regions highlighted", formatter_class=formatter_class)

    parser_cds_group.add_argument("-c", "--df_cosmic", type=str, help="COSMIC mutations() dataframe file path", required=True)
    parser_cds_group.add_argument("-s", "--df_cds", type=str, help="CDS region file path (with columns: gene, CDS, start, end)", required=True)

    parser_cds_group.add_argument("-o", "--out_dir", type=str, help="Output directory for plot", default='../out')

    parser_cds_group.set_defaults(func=co.cds_group)

    # priority_muts: returns the shared sequences library dataframe with priority mutations
    parser_cosmic_priority_muts = subparsers_cosmic.add_parser("priority_muts", help="Identify priority mutations in shared pegRNA library", description="Identify priority mutations in shared pegRNA library", formatter_class=formatter_class)

    parser_cosmic_priority_muts.add_argument("-P", "--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path", required=True)
    parser_cosmic_priority_muts.add_argument("-c", "--df_cosmic", type=str, help="COSMIC mutations() dataframe file path", required=True)

    parser_cosmic_priority_muts.add_argument("-o", "--dir", type=str, help="Output directory", default='../out')
    parser_cosmic_priority_muts.add_argument("-f", "--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_shared_mutations.csv')

    parser_cosmic_priority_muts.set_defaults(func=co.priority_muts)

    # priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    parser_cosmic_priority_edits = subparsers_cosmic.add_parser("priority_edits", help="Identify clinically-relevant prime edits from shared pegRNA sequences", description="Identify clinically-relevant prime edits from shared pegRNA sequences", formatter_class=formatter_class)

    parser_cosmic_priority_edits.add_argument("-p", "--pegRNAs", type=str, help="pegRNAs library dataframe file path", required=True)
    parser_cosmic_priority_edits.add_argument("-P", "--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path", required=True)
    parser_cosmic_priority_edits.add_argument("-c", "--df_cosmic", type=str, help="COSMIC mutations() dataframe file path", required=True)
    
    parser_cosmic_priority_edits.add_argument("-o", "--dir", type=str, help="Output directory", default='../out')
    parser_cosmic_priority_edits.add_argument("-f", "--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_priority.csv')

    parser_cosmic_priority_edits.set_defaults(func=co.priority_edits)

    # editor_mutations(): returns and plots editor accessible COSMIC mutations
    parser_editor_muts = subparsers_cosmic.add_parser("editor_mutations", help="Plot editor-accessible COSMIC mutations using BESCAN library", description="Plot editor-accessible COSMIC mutations using BESCAN library", formatter_class=formatter_class)

    parser_editor_muts.add_argument("-c", "--df_cosmic", type=str, help="COSMIC mutations() dataframe file path", required=True)
    parser_editor_muts.add_argument("-b", "--df_bescan", type=str, help="BESCAN sgRNA library dataframe file path", required=True)

    parser_editor_muts.add_argument("-o", "--out_dir", type=str, help="Output directory for plots or results", default='../out')

    parser_editor_muts.set_defaults(func=co.editor_mutations)

    '''
    edms.data.cvar:
    - mutations(): returns ClinVar mutations dataframe for a given gene
    - priority_muts: returns the shared sequences library dataframe with priority mutations
    - priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    '''
    parser_cvar = subparsers.add_parser("cvar", help="ClinVar Database", description="ClinVar Database", formatter_class=formatter_class)
    subparsers_cvar = parser_cvar.add_subparsers()

    # mutations(): returns ClinVar mutations dataframe for a given gene
    parser_cvar_mutations = subparsers_cvar.add_parser("mutations", help="Extract ClinVar mutations", description="Extract ClinVar mutations", formatter_class=formatter_class)

    parser_cvar_mutations.add_argument("-d", "--df", type=str, help="Input file path", required=True)
    parser_cvar_mutations.add_argument("-g", "--gene_name", type=str, help="Gene name", required=True)

    parser_cvar_mutations.add_argument("-o", "--dir", type=str, help="Output directory", default='../out')
    parser_cvar_mutations.add_argument("-f", "--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_cvar_mutations.csv')

    parser_cvar_mutations.set_defaults(func=cvar.mutations)

    # priority_muts: returns the shared sequences library dataframe with priority mutations
    parser_cvar_priority_muts = subparsers_cvar.add_parser("priority_muts", help="Identify priority mutations in shared pegRNA library", formatter_class=formatter_class)

    parser_cvar_priority_muts.add_argument("-P", "--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path", required=True)
    parser_cvar_priority_muts.add_argument("-c", "--df_clinvar", type=str, help="ClinVar mutations() dataframe file path", required=True)

    parser_cvar_priority_muts.add_argument("-o", "--dir", type=str, help="Output directory", default='../out')
    parser_cvar_priority_muts.add_argument("-f", "--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_shared_mutations.csv')

    parser_cvar_priority_muts.set_defaults(func=cvar.priority_muts)

    # priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    parser_cvar_priority_edits = subparsers_cvar.add_parser("priority_edits", help="Identify clinically-relevant prime edits from shared pegRNA sequences", description="Identify clinically-relevant prime edits from shared pegRNA sequences", formatter_class=formatter_class)

    parser_cvar_priority_edits.add_argument("-p", "--pegRNAs", type=str, help="pegRNAs library dataframe file path", required=True)
    parser_cvar_priority_edits.add_argument("-P", "--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path", required=True)
    parser_cvar_priority_edits.add_argument("-c", "--df_clinvar", type=str, help="ClinVar mutations() dataframe file path", required=True)

    parser_cvar_priority_edits.add_argument("-o", "--dir", type=str, help="Output directory", default='../out')
    parser_cvar_priority_edits.add_argument("-f", "--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_priority.csv')

    parser_cvar_priority_edits.set_defaults(func=cvar.priority_edits)

    '''
    edms.data.cBioPortal:
    - mutations(): Process cBioPortal mutation data for a single protein-coding gene.
    '''
    parser_cBP = subparsers.add_parser("cBioPortal", help="cBioPortal Database", description="cBioPortal Database", formatter_class=formatter_class)
    subparsers_cBP = parser_cBP.add_subparsers()

    # mutations(): returns cBioPortal mutations dataframe for a given gene
    parser_cBP_mutations = subparsers_cBP.add_parser("mutations", help="Process cBioPortal mutation data for a single protein-coding gene. Retrieve cBioPortal mutation data from https://www.cbioportal.org/ (GENIE Cohort v18.0-public).", description="Process cBioPortal mutation data for a single protein-coding gene. Retrieve cBioPortal mutation data from https://www.cbioportal.org/ (GENIE Cohort v18.0-public).", formatter_class=formatter_class)
    
    parser_cBP_mutations.add_argument("-d", "--df", type=str, help="Path to cBioPortal mutation data TSV file", required=True)
    parser_cBP_mutations.add_argument("-w", "--wt", type=str, help="Wild-type amino acid sequence", required=True)
    parser_cBP_mutations.add_argument("-C", "--config", type=bool, help="Save to configuration directory", default=True)
    parser_cBP_mutations.add_argument("-o", "--dir", type=str, help="Output directory", default='../out')
    parser_cBP_mutations.add_argument("-f", "--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_cBioPortal_mutations.csv')

    parser_cBP_mutations.set_defaults(func=cBP.mutations)

    '''
    edms.data.uniprot:
    - retrieve(): Download a UniProt flat file entry via the REST API and save it to config and/or specified directory.
    '''
    parser_uniprot = subparsers.add_parser("uniprot", help="UniProt Database", description="UniProt Database", formatter_class=formatter_class)
    subparsers_uniprot = parser_uniprot.add_subparsers()

    # retrieve(): Download a UniProt flat file entry via the REST API and save it to config and/or specified directory.
    parser_uniprot_retrieve = subparsers_uniprot.add_parser("retrieve", help="Download a UniProt flat file entry via the REST API and save it to config and/or specified directory.", description="Download a UniProt flat file entry via the REST API and save it to config and/or specified directory.", formatter_class=formatter_class)
    
    parser_uniprot_retrieve.add_argument("-a", "--accession", type=str, help="UniProt accession (e.g. P55317)", required=True)
    parser_uniprot_retrieve.add_argument("-o", "--dir", type=str, help="Save flat file to specified directory (Default: None -> only save to config directory)", default=argparse.SUPPRESS)
    parser_uniprot_retrieve.add_argument("-n", "--no_config", action="store_false", dest='config', help="Do not save to configuration directory", default=True)
    parser_uniprot_retrieve.add_argument("-b", "--base_url", type=str, help="Base URL for UniProt REST API (Default: https://rest.uniprot.org/uniprotkb)", default="https://rest.uniprot.org/uniprotkb")


    parser_uniprot_retrieve.set_defaults(func=uniprot.retrieve)

    # draw_ss_track(): Draw a secondary structure track
    parser_uniprot_draw_ss_track = subparsers_uniprot.add_parser("ss", help="Draw a secondary structure track", description="Draw a secondary structure track", formatter_class=formatter_class)

    parser_uniprot_draw_ss_track.add_argument("-i", "--df", type=str, help="UniProt accession (if saved to ~/.config/edms/UniProt) or file path for UniProt flat file. See edms.dat.uniprot.retrieve() or edms uniprot retrieve -h for more information.", required=True)
    
    parser_uniprot_draw_ss_track.add_argument("-y", "--y", type=float, help="Vertical position of the track", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-ts", "--track_start", type=float, help="Optional start bound for the track; if None, determined from data", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-te", "--track_end", type=float, help="Optional end bound for the track; if None, determined from data", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-hh", "--helix_height", type=float, help="Height for helices; if None, defaults to `height`", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-sh", "--strand_height", type=float, help="Height for strands; if None, defaults to `height`", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-lh", "--loop_height", type=float, help="Height for loops; if None, defaults to `height`", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-xp", "--xpad", type=int, help="Horizontal padding on the left and right of the track", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-hc", "--helix_color", type=str, help="Color for helices", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-sc", "--strand_color", type=str, help="Color for strands", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-lc", "--loop_color", type=str, help="Color for loops", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-ec", "--edgecolor", type=str, help="Color of the edges of helices and strands", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-hlw", "--helix_lw", type=float, help="Line width for helices", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-slw", "--strand_lw", type=float, help="Line width for strands", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-llw", "--loop_lw", type=float, help="Line width for loops", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-lml", "--loop_min_len", type=float, help="Minimum length of loops to draw (for auto-drawn loops between elements)", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-z", "--base_zorder", type=int, help="Z-order for the first element; subsequent elements decrease by zorder_step to create a stacking effect", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-hlp", "--helix_loop_pad", type=float, help="Padding to apply to the start/end of helices when auto-drawing loops between elements, to make the sine wave connect more nicely; applied as a subtraction from the end of the current element and an addition to the start of the next element when determining loop start and end positions", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-slp", "--strand_loop_pad", type=float, help="Padding to apply to the start/end of strands when auto-drawing loops between elements, to make the sine wave connect more nicely; applied as a subtraction from the end of the current element and an addition to the start of the next element when determining loop start and end positions", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-fs","--figsize", type=parse_tuple_float, help="Size of the figure to create (width, height)", default=(10, 2))
    parser_uniprot_draw_ss_track.add_argument("-o", "--dir", type=str, help="Output directory", default='../out')
    parser_uniprot_draw_ss_track.add_argument("-f", "--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_uniprot_ss_track.all')
    parser_uniprot_draw_ss_track.add_argument("-d", "--dpi", type=int, help="Figure dpi (Default: 600)", default=argparse.SUPPRESS)
    parser_uniprot_draw_ss_track.add_argument("-s", "--show", action="store_true", help="Show plot", default=False)

    parser_uniprot_draw_ss_track.set_defaults(func=uniprot.draw_ss_track)

    '''
    edms.data.pdb:
    - retrieve(): Download a PDB or mmCIF file from the RCSB PDB REST API.
    '''
    parser_pdb = subparsers.add_parser("pdb", help="PDB Database", description="PDB Database", formatter_class=formatter_class)
    subparsers_pdb = parser_pdb.add_subparsers()

    # retrieve(): Download a PDB or PDBx/mmCIF file from the RCSB PDB REST API.
    parser_pdb_retrieve = subparsers_pdb.add_parser("retrieve", help="Download a PDB or PDBx/mmCIF file from the RCSB PDB REST API.", description="Download a PDB or PDBx/mmCIF file from the RCSB PDB REST API.", formatter_class=formatter_class)
    
    parser_pdb_retrieve.add_argument("-i", "--id", type=str, help="4-character PDB accession (e.g. 8VG1, case-insensitive)", required=True)
    parser_pdb_retrieve.add_argument("-s", "--suf", type=str, help="Output filename suffix (e.g., .pdb or .cif (PDBx/mmCIF))", default=argparse.SUPPRESS)
    parser_pdb_retrieve.add_argument("-o", "--dir", type=str, help="Save file to specified directory (Default: None -> only save to config directory)", default=argparse.SUPPRESS)
    parser_pdb_retrieve.add_argument("-n", "--no_config", action="store_false", dest='config', help="Do not save to configuration directory", default=True)
    parser_pdb_retrieve.add_argument("-b", "--base_url", type=str, help="Base URL for RCSB PDB REST API (Default: https://files.rcsb.org/download/)", default="https://files.rcsb.org/download/")

    parser_pdb_retrieve.set_defaults(func=pdb.retrieve)

    '''
    edms.data.dssp:
    - retrieve(): Retrieve DSSP data for a given PDB ID.
    '''
    parser_dssp = subparsers.add_parser("dssp", help="Define Secondary Structure of Proteins", description="Define Secondary Structure of Proteins", formatter_class=formatter_class)

    subparsers_dssp = parser_dssp.add_subparsers()
    
    # retrieve(): Retrieve DSSP data for a given PDB ID.
    parser_dssp_retrieve = subparsers_dssp.add_parser("retrieve", help="Retrieve DSSP data for a given PDB ID.", description="Retrieve DSSP data for a given PDB ID.", formatter_class=formatter_class)

    parser_dssp_retrieve.add_argument("-i", "--id_or_file", type=str, help="4-character PDB accession (e.g. 8VG1, case-insensitive) or file path for PDB structure file. See edms.data.dssp.retrieve() for more information.", required=True)
    parser_dssp_retrieve.add_argument("-o", "--dir", type=str, help="Save file to specified directory; otherwise, save to configuration directory (~/.config/edms/DSSP/).", default=argparse.SUPPRESS)

    parser_dssp_retrieve.set_defaults(func=dssp.retrieve)

    # plot_ss_color_key(): Plot a color key for DSSP secondary structure assignments.
    parser_dssp_plot_ss_color_key = subparsers_dssp.add_parser("color_key", help="Plot a color key for DSSP secondary structure assignments.", description="Plot a color key for DSSP secondary structure assignments.", formatter_class=formatter_class)

    parser_dssp_plot_ss_color_key.add_argument("-n", "--ncols", type=int, help="Number of columns in the key", default=argparse.SUPPRESS)
    parser_dssp_plot_ss_color_key.add_argument("-fs", "--figsize", type=parse_tuple_float, help="Figure size (width, height)", default=argparse.SUPPRESS)
    parser_dssp_plot_ss_color_key.add_argument("-t", "--title", type=str, help="Figure title", default=argparse.SUPPRESS)
    parser_dssp_plot_ss_color_key.add_argument("-ts", "--title_size", type=int, help="Title font size", default=argparse.SUPPRESS)
    parser_dssp_plot_ss_color_key.add_argument("-tw", "--title_weight", type=str, help="Title font weight (e.g., 'bold', 'normal')", default=argparse.SUPPRESS)
    parser_dssp_plot_ss_color_key.add_argument("-F", "--fontsize", type=int, help="Label font size", default=argparse.SUPPRESS)
    parser_dssp_plot_ss_color_key.add_argument("-sw", "--swatch_width", type=float, help="Width of color box", default=argparse.SUPPRESS)
    parser_dssp_plot_ss_color_key.add_argument("-sh", "--swatch_height", type=float, help="Height of color box", default=argparse.SUPPRESS)
    parser_dssp_plot_ss_color_key.add_argument("-a", "--alpha", type=float, help="Transparency of color boxes", default=argparse.SUPPRESS)
    parser_dssp_plot_ss_color_key.add_argument("-o", "--dir", type=str, help="Directory to save the figure; if not provided, the figure is not saved", default='.')
    parser_dssp_plot_ss_color_key.add_argument("-f", "--file", type=str, help="Filename to save the figure; if not provided, the figure is not saved", default='dssp_color_key.all')
    parser_dssp_plot_ss_color_key.add_argument("-d", "--dpi", type=int, help="Resolution for saving the figure (Default: 1200)", default=argparse.SUPPRESS)
    parser_dssp_plot_ss_color_key.add_argument("-nt", "--not_transparent", dest='transparent', action="store_false", help="Don't save plot with transparent background", default=True)
    parser_dssp_plot_ss_color_key.add_argument("-s", "--show", action="store_true", help="Show plot", default=False)

    parser_dssp_plot_ss_color_key.set_defaults(func=dssp.plot_ss_color_key)