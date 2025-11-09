'''
src/edms/dat/cli.py             Command Line Interface for EDMS Database module
├── __init__.py                 Initializer
├── cosmic.py                   COSMIC module
├── cvar.py                     ClinVar module
└── ncbi.py                     NCBI module
'''
import argparse
import datetime
import sys # might use later
from rich import print as rprint # might use later

from . import cosmic as co, cvar
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
    edms.dat.cosmic:
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

    parser_cosmic_mutations.add_argument("--df", type=str, help="Input file path", required=True)

    parser_cosmic_mutations.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cosmic_mutations.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_cosmic_mutations.csv')

    parser_cosmic_mutations.set_defaults(func=co.mutations)
    
    # cds_group(): plot COSMIC mutations histogram with CDS regions highlighted in different colors
    parser_cds_group = subparsers_cosmic.add_parser("cds_group", help="Plot COSMIC mutation histogram with CDS regions highlighted", description="Plot COSMIC mutation histogram with CDS regions highlighted", formatter_class=formatter_class)

    parser_cds_group.add_argument("--df_cosmic", type=str, help="COSMIC mutations() file path", required=True)
    parser_cds_group.add_argument("--df_cds", type=str, help="CDS region file path (with columns: gene, CDS, start, end)", required=True)

    parser_cds_group.add_argument("--out_dir", type=str, help="Output directory for plot",default='../out')

    parser_cds_group.set_defaults(func=co.cds_group)

    # priority_muts: returns the shared sequences library dataframe with priority mutations
    parser_cosmic_priority_muts = subparsers_cosmic.add_parser("priority_muts", help="Identify priority mutations in shared pegRNA library", description="Identify priority mutations in shared pegRNA library", formatter_class=formatter_class)

    parser_cosmic_priority_muts.add_argument("--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path", required=True)
    parser_cosmic_priority_muts.add_argument("--df_cosmic", type=str, help="COSMIC mutations() dataframe file path",required=True)

    parser_cosmic_priority_muts.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cosmic_priority_muts.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_shared_mutations.csv')

    parser_cosmic_priority_muts.set_defaults(func=co.priority_muts)

    # priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    parser_cosmic_priority_edits = subparsers_cosmic.add_parser("priority_edits", help="Identify clinically-relevant prime edits from shared pegRNA sequences", description="Identify clinically-relevant prime edits from shared pegRNA sequences", formatter_class=formatter_class)

    parser_cosmic_priority_edits.add_argument("--pegRNAs", type=str, help="pegRNAs library dataframe file path",required=True)
    parser_cosmic_priority_edits.add_argument("--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path",required=True)
    parser_cosmic_priority_edits.add_argument("--df_cosmic", type=str, help="COSMIC mutations() dataframe file path",required=True)
    
    parser_cosmic_priority_edits.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cosmic_priority_edits.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_priority.csv')

    parser_cosmic_priority_edits.set_defaults(func=co.priority_edits)

    # editor_mutations(): returns and plots editor accessible COSMIC mutations
    parser_editor_muts = subparsers_cosmic.add_parser("editor_mutations", help="Plot editor-accessible COSMIC mutations using BESCAN library", description="Plot editor-accessible COSMIC mutations using BESCAN library", formatter_class=formatter_class)

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
    parser_cvar = subparsers.add_parser("cvar", help="ClinVar Database", description="ClinVar Database", formatter_class=formatter_class)
    subparsers_cvar = parser_cvar.add_subparsers()

    # mutations(): returns ClinVar mutations dataframe for a given gene
    parser_cvar_mutations = subparsers_cvar.add_parser("mutations", help="Extract ClinVar mutations", description="Extract ClinVar mutations", formatter_class=formatter_class)

    parser_cvar_mutations.add_argument("--df", type=str, help="Input file path", required=True)
    parser_cvar_mutations.add_argument("--gene_name", type=str, help="Gene name", required=True)

    parser_cvar_mutations.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cvar_mutations.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_cvar_mutations.csv')

    parser_cvar_mutations.set_defaults(func=cvar.mutations)

    # priority_muts: returns the shared sequences library dataframe with priority mutations
    parser_cvar_priority_muts = subparsers_cvar.add_parser("priority_muts", help="Identify priority mutations in shared pegRNA library", formatter_class=formatter_class)

    parser_cvar_priority_muts.add_argument("--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path", required=True)
    parser_cvar_priority_muts.add_argument("--df_clinvar", type=str, help="ClinVar mutations() dataframe file path",required=True)

    parser_cvar_priority_muts.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cvar_priority_muts.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_shared_mutations.csv')

    parser_cvar_priority_muts.set_defaults(func=cvar.priority_muts)

    # priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    parser_cvar_priority_edits = subparsers_cvar.add_parser("priority_edits", help="Identify clinically-relevant prime edits from shared pegRNA sequences", description="Identify clinically-relevant prime edits from shared pegRNA sequences", formatter_class=formatter_class)

    parser_cvar_priority_edits.add_argument("--pegRNAs", type=str, help="pegRNAs library dataframe file path",required=True)
    parser_cvar_priority_edits.add_argument("--pegRNAs_shared", type=str, help="Shared pegRNAs library dataframe file path",required=True)
    parser_cvar_priority_edits.add_argument("--df_clinvar", type=str, help="ClinVar mutations() dataframe file path",required=True)

    parser_cvar_priority_edits.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_cvar_priority_edits.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pegRNAs_priority.csv')

    parser_cvar_priority_edits.set_defaults(func=cvar.priority_edits)