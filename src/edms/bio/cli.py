''' src/edms/bio/cli.py         Command Line Interface for EDMS Biology module
├── clone.py                    Clone module
├── fastq.py                    FASTQ module
├── ngs.py                      Next-generation Sequencing module
├── pe.py                       Prime Editing module
├── pegLIT.py                   pegRNA Linker module
├── primedesign.py              PrimeDesign module
├── qPCR.py                     qPCR module
├── sanger.py                   Sanger Sequencing module
├── signature.py                Signature module
└── transfect.py                Transfection module
'''
import argparse
import datetime
import sys
from rich import print as rprint

from . import ngs, sanger, clone as cl, fastq as fq, pe, qPCR, transfect as tf
from ..utils import parse_tuple_int, parse_tuple_float
from ..gen.cli import add_common_plot_cat_args, add_common_plot_heat_args, add_common_plot_scat_args, add_common_plot_stack_args, add_common_plot_vol_args

def add_subparser(subparsers, formatter_class=None):
    """
    add_subparser(): Attach all bio-related subparsers to the top-level CLI.

    Parameters:
    subparsers (argparse._SubParsersAction): The subparsers object to attach the bio subparsers to.
    formatter_class (type, optional): The formatter class to use for the subparsers.
    """
    if formatter_class is None:
        # fall back to basic formatter to avoid circular imports
        formatter_class = argparse.HelpFormatter

    '''
    edms.bio.ngs:
    - pcrs(): generates NGS PCR plan automatically
    - umis(): determine the ug of gDNA and reads required for genotyping with UMIs.
    - compute_distance_matrix(): compute pairwise Hamming distance matrix for a list of sequences stored in a dataframe
    '''
    parser_ngs = subparsers.add_parser("ngs", help="Next generation sequencing", description="Next generation sequencing", formatter_class=formatter_class)
    subparsers_ngs = parser_ngs.add_subparsers()

    # pcrs(): generates NGS PCR plan automatically
    parser_ngs_pcrs = subparsers_ngs.add_parser("pcrs", help="Plan NGS PCRs", description="Plan NGS PCRs", formatter_class=formatter_class)

    # pcrs(): Core parameters
    parser_ngs_pcrs.add_argument("--df", help="Input file", type=str, required=True)
    
    parser_ngs_pcrs.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_ngs_pcrs.add_argument("--file", help="Output file name (.xlsx)", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_NGS_plan.xlsx')
    parser_ngs_pcrs.add_argument("--ultra", help="Using NEB Ultra II reagents", action="store_true")
    parser_ngs_pcrs.add_argument("--pcr1_cycles", help="Number of cycles for PCR1", type=str, default='30')
    parser_ngs_pcrs.add_argument("--pcr2_cycles", help="Number of cycles for PCR2", type=str, default='8')
    parser_ngs_pcrs.add_argument("--umi_cycles", help="Number of cycles for PCR1", type=str, default='3')
    parser_ngs_pcrs.add_argument("--pcr1.5_Tm",dest='pcr1_5_Tm', help="Annealing temperature for PCR1.5", type=str, default='65')
    parser_ngs_pcrs.add_argument('--pcr1_total_uL', type=int, default=20, help='PCR1 Total reaction volume (uL)')
    parser_ngs_pcrs.add_argument('--pcr2_total_uL', type=int, default=20, help='PCR2 Total reaction volume (uL)')
    parser_ngs_pcrs.add_argument('--mm_x', type=float, default=1.1, help='Master mix multiplier')
    
    # pcrs(): Column names
    parser_ngs_pcrs.add_argument('--gDNA_id_col', default='ID', help='gDNA ID column name')
    parser_ngs_pcrs.add_argument('--pcr1_id_col', default='PCR1 ID', help='PCR1 ID column name')
    parser_ngs_pcrs.add_argument('--pcr1_fwd_col', default='PCR1 FWD', help='PCR1 FWD column name')
    parser_ngs_pcrs.add_argument('--pcr1_rev_col', default='PCR1 REV', help='PCR1 REV column name')
    parser_ngs_pcrs.add_argument('--pcr2_id_col', default='PCR2 ID', help='PCR2 ID column name')
    parser_ngs_pcrs.add_argument('--pcr2_fwd_col', default='PCR2 FWD', help='PCR2 FWD column name')
    parser_ngs_pcrs.add_argument('--pcr2_rev_col', default='PCR2 REV', help='PCR2 REV column name')
    parser_ngs_pcrs.add_argument('--umi_col', default='UMI', help='UMI column name')

    # pcrs(): Stock concentrations
    parser_ngs_pcrs.add_argument('--Q5_mm_x_stock', type=float, default=5, help='Q5 reaction master mix stock (X)')
    parser_ngs_pcrs.add_argument('--dNTP_mM_stock', type=float, default=10, help='dNTP stock concentration (mM)')
    parser_ngs_pcrs.add_argument('--fwd_uM_stock', type=float, default=10, help='Forward primer stock concentration (uM)')
    parser_ngs_pcrs.add_argument('--rev_uM_stock', type=float, default=10, help='Reverse primer stock concentration (uM)')
    parser_ngs_pcrs.add_argument('--Q5_U_uL_stock', type=float, default=2, help='Q5 Polymerase stock (U/uL)')

    # pcrs(): Desired concentrations
    parser_ngs_pcrs.add_argument('--Q5_mm_x_desired', type=float, default=1, help='Q5 reaction master mix desired (X)')
    parser_ngs_pcrs.add_argument('--dNTP_mM_desired', type=float, default=0.2, help='dNTP desired concentration (mM)')
    parser_ngs_pcrs.add_argument('--fwd_uM_desired', type=float, default=0.5, help='Forward primer desired concentration (uM)')
    parser_ngs_pcrs.add_argument('--rev_uM_desired', type=float, default=0.5, help='Reverse primer desired concentration (uM)')
    parser_ngs_pcrs.add_argument('--Q5_U_uL_desired', type=float, default=0.02, help='Q5 Polymerase desired amount (U/uL)')

    parser_ngs_pcrs.set_defaults(func=ngs.pcrs)
    
    # umis(): determine the ug of gDNA and reads required for genotyping with UMIs.
    parser_ngs_umis = subparsers_ngs.add_parser("umis", help="Determine ug of gDNA and reads required for genotyping with UMIs", description="Determine ug of gDNA and reads required for genotyping with UMIs", formatter_class=formatter_class)

    parser_ngs_umis.add_argument("--genotypes", help="# of intended genotypes per sample", type=int, required=True)

    parser_ngs_umis.add_argument("--samples", help="# of samples to be processed (Default: 1", type=int, default=1)
    parser_ngs_umis.add_argument("--cell_coverage", help="Desired coverage per genotype (Default: 1000)", type=int, default=1000)
    parser_ngs_umis.add_argument("--ug_gDNA_per_cell", help="Amount of genomic DNA per cell in micrograms (Default: 6x10^(-6) ug/cell)", type=float, default=6*10**-6)
    parser_ngs_umis.add_argument("--ploidy_per_cell", help="Ploidy level of the cells (Default: 2 for diploid)", type=int, default=2)
    parser_ngs_umis.add_argument("--umi_coverage", help="Average # of reads per UMI (Default: 5)", type=int, default=5)

    parser_ngs_umis.set_defaults(func=ngs.umis)
    
    # hamming_distance_matrix(): compute pairwise Hamming distance matrix for a list of sequences stored in a dataframe
    parser_ngs_hamming = subparsers_ngs.add_parser("hamming", help="Compute pairwise Hamming distance matrix", description="Compute pairwise Hamming distance matrix", formatter_class=formatter_class)
    
    parser_ngs_hamming.add_argument("--df", help="Input file", type=str, required=True)
    parser_ngs_hamming.add_argument("--id", help="ID column name", type=str, required=True)
    parser_ngs_hamming.add_argument("--seqs", help="Sequences column name", type=str, required=True)
    
    parser_ngs_hamming.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_ngs_hamming.add_argument("--file", help="Output file name", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_hamming.csv')
    
    parser_ngs_hamming.set_defaults(func=ngs.hamming_distance_matrix)

    '''
    edms.bio.sanger:
    - pcrs(): generates Sanger PCR plan automatically
    '''
    parser_sanger = subparsers.add_parser("sanger", help="Sanger sequencing", description="Sanger sequencing", formatter_class=formatter_class)
    subparsers_sanger = parser_sanger.add_subparsers()

    # pcrs(): generates Sanger PCR plan automatically
    parser_sanger_pcrs = subparsers_sanger.add_parser("pcrs", help="Plan Sanger PCRs", description="Plan Sanger PCRs", formatter_class=formatter_class)

    # pcrs(): Core parameters
    parser_sanger_pcrs.add_argument("--df", help="Input file", type=str, required=True)
    parser_sanger_pcrs.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_sanger_pcrs.add_argument("--file", help="Output file name (.xlsx)", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_Sanger_plan.xlsx')
    parser_sanger_pcrs.add_argument("--cycles", help="Number of cycles for PCR1", type=str, default='30')
    parser_sanger_pcrs.add_argument("--ultra", help="Using NEB Ultra II reagents", action="store_true")
    parser_sanger_pcrs.add_argument('--total_uL', type=int, default=20, help='Total reaction volume (uL)')
    parser_sanger_pcrs.add_argument('--mm_x', type=float, default=1.1, help='Master mix multiplier')
    
    # pcrs(): Column names
    parser_sanger_pcrs.add_argument('--gDNA_id_col', default='ID', help='gDNA ID column name')
    parser_sanger_pcrs.add_argument('--pcr1_id_col', default='PCR1 ID', help='PCR1 ID column name')
    parser_sanger_pcrs.add_argument('--pcr1_fwd_col', default='PCR1 FWD', help='PCR1 FWD column name')
    parser_sanger_pcrs.add_argument('--pcr1_rev_col', default='PCR1 REV', help='PCR1 REV column name')
   
    # pcrs(): Stock concentrations
    parser_sanger_pcrs.add_argument('--Q5_mm_x_stock', type=float, default=5, help='Q5 reaction master mix stock (X)')
    parser_sanger_pcrs.add_argument('--dNTP_mM_stock', type=float, default=10, help='dNTP stock concentration (mM)')
    parser_sanger_pcrs.add_argument('--fwd_uM_stock', type=float, default=10, help='Forward primer stock concentration (uM)')
    parser_sanger_pcrs.add_argument('--rev_uM_stock', type=float, default=10, help='Reverse primer stock concentration (uM)')
    parser_sanger_pcrs.add_argument('--Q5_U_uL_stock', type=float, default=2, help='Q5 Polymerase stock (U/uL)')

    # pcrs(): Desired concentrations
    parser_sanger_pcrs.add_argument('--Q5_mm_x_desired', type=float, default=1, help='Q5 reaction master mix desired (X)')
    parser_sanger_pcrs.add_argument('--dNTP_mM_desired', type=float, default=0.2, help='dNTP desired concentration (mM)')
    parser_sanger_pcrs.add_argument('--fwd_uM_desired', type=float, default=0.5, help='Forward primer desired concentration (uM)')
    parser_sanger_pcrs.add_argument('--rev_uM_desired', type=float, default=0.5, help='Reverse primer desired concentration (uM)')
    parser_sanger_pcrs.add_argument('--Q5_U_uL_desired', type=float, default=0.02, help='Q5 Polymerase desired amount (U/uL)')

    parser_sanger_pcrs.set_defaults(func=sanger.pcrs)
    
    '''
    edms.bio.clone
    - sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
    - epegRNAs(): design GG cloning oligonucleotides for prime editing epegRNAs
    - ngRNAs(): design GG cloning oligonucleotides for prime editing ngRNAs
    - epegRNA_pool(): makes twist oligonucleotides for prime editing
    - pcr_sim(): returns dataframe with simulated pcr product
    - off_targets(): Find off-target sequences for a list of sequences using pairwise alignment.
    '''
    parser_clone = subparsers.add_parser("clone", help="Molecular cloning", description="Molecular cloning", formatter_class=formatter_class)
    subparsers_clone = parser_clone.add_subparsers()

    # sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
    parser_clone_sgRNAs = subparsers_clone.add_parser("sgRNAs", help="Design GG oligos for sgRNAs (cutting or BE)", description="Design GG oligos for sgRNAs (cutting or BE)", formatter_class=formatter_class)

    parser_clone_sgRNAs.add_argument("--df", type=str, help="Input file path",required=True)
    parser_clone_sgRNAs.add_argument("--id", type=str, help="Column name for unique sgRNA identifier",required=True)

    parser_clone_sgRNAs.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_clone_sgRNAs.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_sgRNAs.csv')

    parser_clone_sgRNAs.add_argument("--dont_tG", dest="tG", default=True, action="store_false", help="Don't add 5' G to spacer if needed")
    parser_clone_sgRNAs.add_argument("--dont_order", dest="order", default=True, action="store_false", help="Don't format output for ordering oligos")
    parser_clone_sgRNAs.add_argument("--spacer", type=str, default="Spacer_sequence", help="Column name for spacer sequence")
    parser_clone_sgRNAs.add_argument("--t5", type=str, default="CACC", help="Top oligo 5' overhang")
    parser_clone_sgRNAs.add_argument("--t3", type=str, default="", help="Top oligo 3' overhang")
    parser_clone_sgRNAs.add_argument("--b5", type=str, default="AAAC", help="Bottom oligo 5' overhang (revcom)")
    parser_clone_sgRNAs.add_argument("--b3", type=str, default="", help="Bottom oligo 3' overhang (revcom)")

    parser_clone_sgRNAs.set_defaults(func=cl.sgRNAs)

    # epegRNAs(): design GG cloning oligonucleotides for prime editing epegRNAs
    parser_clone_epegRNAs = subparsers_clone.add_parser("epegRNAs", help="Design GG oligos for epegRNAs", description="Design GG oligos for epegRNAs", formatter_class=formatter_class)

    parser_clone_epegRNAs.add_argument("--df", type=str, help="Input file path", required=True)
    parser_clone_epegRNAs.add_argument("--id", type=str, help="Column name for unique sequence identifier",required=True)

    parser_clone_epegRNAs.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_clone_epegRNAs.add_argument("--file", help="Output file name", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_epegRNAs.csv')

    parser_clone_epegRNAs.add_argument("--dont_tG", dest="tG", default=True, action="store_false", help="Don't add 5' G to spacer if needed")
    parser_clone_epegRNAs.add_argument("--dont_order", dest="order", default=True, action="store_false", help="Don't format output for ordering oligos")
    parser_clone_epegRNAs.add_argument("--order_scaffold", action="store_true", help="Include scaffold sequence in the oligo order")
    parser_clone_epegRNAs.add_argument("--dont_make_extension", dest="make_extension", default=True, action="store_false", help="Don't build extension from RTT, PBS, and linker")
    parser_clone_epegRNAs.add_argument("--spacer", type=str, default="Spacer_sequence", help="Column name for spacer sequence")
    parser_clone_epegRNAs.add_argument("--spacer_t5", type=str, default="CACC", help="Top 5' overhang for spacer")
    parser_clone_epegRNAs.add_argument("--spacer_t3", type=str, default="GTTTAAGAGC", help="Top 3' overhang for spacer")
    parser_clone_epegRNAs.add_argument("--spacer_b5", type=str, default="", help="Bottom 5' overhang for spacer")
    parser_clone_epegRNAs.add_argument("--spacer_b3", type=str, default="", help="Bottom 3' overhang for spacer")
    parser_clone_epegRNAs.add_argument("--extension", type=str, default="Extension_sequence", help="Column name for extension sequence")
    parser_clone_epegRNAs.add_argument("--extension_t5", type=str, default="", help="Top 5' overhang for extension")
    parser_clone_epegRNAs.add_argument("--extension_t3", type=str, default="", help="Top 3' overhang for extension")
    parser_clone_epegRNAs.add_argument("--extension_b5", type=str, default="CGCG", help="Bottom 5' overhang for extension")
    parser_clone_epegRNAs.add_argument("--extension_b3", type=str, default="GCACCGACTC", help="Bottom 3' overhang for extension")
    parser_clone_epegRNAs.add_argument("--RTT", type=str, default="RTT_sequence", help="Column name for RTT (reverse transcriptase template)")
    parser_clone_epegRNAs.add_argument("--PBS", type=str, default="PBS_sequence", help="Column name for PBS (primer binding site)")
    parser_clone_epegRNAs.add_argument("--linker", type=str, default="Linker_sequence", help="Column name for linker")

    parser_clone_epegRNAs.set_defaults(func=cl.epegRNAs)
    
    # ngRNAs(): design GG cloning oligonucleotides for prime editing ngRNAs
    parser_clone_ngRNAs = subparsers_clone.add_parser("ngRNAs", help="Design GG oligos for ngRNAs", description="Design GG oligos for ngRNAs", formatter_class=formatter_class)

    parser_clone_ngRNAs.add_argument("--df", type=str, help="Input file path", required=True)
    parser_clone_ngRNAs.add_argument("--id", type=str, help="Column name for unique sequence identifier",required=True)

    parser_clone_ngRNAs.add_argument("--dir", help="Output directory path", type=str, default='../out')
    parser_clone_ngRNAs.add_argument("--file", help="Output file name", type=str, default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_epegRNAs.csv')
    
    parser_clone_ngRNAs.add_argument("--dont_tG", dest="tG", default=True, action="store_false", help="Don't add 5' G to spacer if needed")
    parser_clone_ngRNAs.add_argument("--dont_order", dest="order", default=True, action="store_false", help="Don't format output for ordering oligos")
    parser_clone_ngRNAs.add_argument("--order_scaffold", action="store_true", help="Include scaffold sequence in the oligo order")
    parser_clone_ngRNAs.add_argument("--spacer", type=str, default="Spacer_sequence", help="Column name for spacer sequence")
    parser_clone_ngRNAs.add_argument("--spacer_t5", type=str, default="CACC", help="Top strand 5' overhang")
    parser_clone_ngRNAs.add_argument("--spacer_t3", type=str, default="GTTTAAGAGC", help="Top strand 3' overhang")
    parser_clone_ngRNAs.add_argument("--spacer_b5", type=str, default="", help="Bottom strand 5' overhang")
    parser_clone_ngRNAs.add_argument("--spacer_b3", type=str, default="", help="Bottom strand 3' overhang")

    parser_clone_ngRNAs.set_defaults(func=cl.ngRNAs)

    # epegRNA_pool(): makes twist oligonucleotides for prime editing
    parser_clone_epegRNA_pool = subparsers_clone.add_parser("epegRNA_pool", help="Design GG oligos for pooled epegRNAs", description="Design GG oligos for pooled epegRNAs", formatter_class=formatter_class)
    
    parser_clone_epegRNA_pool.add_argument("--df", type=str, help="Input file path", required=True)

    parser_clone_epegRNA_pool.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_clone_epegRNA_pool.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_epegRNA_pool.csv')

    parser_clone_epegRNA_pool.add_argument("--tG", action="store_true", help="Add 5' G to spacers if needed")
    parser_clone_epegRNA_pool.add_argument("--make_extension", action="store_true", help="Build extension from RTT, PBS, and linker")
    parser_clone_epegRNA_pool.add_argument("--UMI_df", type=str, help="UMI sequences file path", default=argparse.SUPPRESS)
    parser_clone_epegRNA_pool.add_argument("--PCR_df", type=str, help="PCR primer and subpool barcode file path", default=argparse.SUPPRESS)
    parser_clone_epegRNA_pool.add_argument("--RE_type_IIS_df", type=str, help="RE Type IIS file path", default=argparse.SUPPRESS)
    parser_clone_epegRNA_pool.add_argument("--UMI_i", type=int, help="UMI start index (Default: 0)", default=0)
    parser_clone_epegRNA_pool.add_argument("--enzymes", type=str, nargs="+", help="List of Type IIS restriction enzymes to check for (Default: Esp3I)", default='Esp3I')
    parser_clone_epegRNA_pool.add_argument("--barcode", type=str, help="subpool barcode column name (Default: Barcode)", default='Barcode')
    parser_clone_epegRNA_pool.add_argument("--barcode_i", type=int, help="subpool barcode start index (Default: 0)", default=0)
    parser_clone_epegRNA_pool.add_argument("--fwd_barcode_t5", type=str, default="Forward Barcode", help="Forward barcode column name")
    parser_clone_epegRNA_pool.add_argument("--rev_barcode_t3", type=str, default="Reverse Barcode", help="Reverse barcode column name")
    parser_clone_epegRNA_pool.add_argument("--Esp3I_hU6", type=str, default="Esp3I_hU6", help="Esp3I_hU6 column name")
    parser_clone_epegRNA_pool.add_argument("--tevopreQ1_Esp3I", type=str, default="tevopreQ1_Esp3I", help="tevopreQ1_Esp3I column name")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_spacer", type=str, default="Spacer_sequence", help="epegRNA spacer column")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_scaffold", type=str, default="Scaffold_sequence", help="epegRNA scaffold column")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_extension", type=str, default="Extension_sequence", help="epegRNA extension column")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_RTT", type=str, default="RTT_sequence", help="epegRNA RTT column name")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_PBS", type=str, default="PBS_sequence", help="epegRNA PBS column name")
    parser_clone_epegRNA_pool.add_argument("--epegRNA_linker", type=str, default="Linker_sequence", help="epegRNA Linker column name")

    parser_clone_epegRNA_pool.set_defaults(func=cl.epegRNA_pool) 
    
    # umi(): generates unique molecular identifiers (UMIs) of specified length, GC content, and Hamming distance
    parser_clone_umi = subparsers_clone.add_parser("umi", help="Generate unique molecular identifiers (UMIs) of specified length, GC content, and Hamming distance", description="Generate unique molecular identifiers (UMIs) of specified length, GC content, and Hamming distance", formatter_class=formatter_class)

    parser_clone_umi.add_argument("--length", type=int, help="Length of UMI (Default: 15)", default=15)
    parser_clone_umi.add_argument("--GC_fract", type=parse_tuple_float, default=(0.4,0.6), help="Pair of GC content boundaries written as fractions (Default: 0.4,0.6)")
    parser_clone_umi.add_argument("--hamming", type=int, help="Minimum Hamming distance between UMIs (Default: 4)", default=4)
    parser_clone_umi.add_argument("--nrows", type=int, help="# of UMIs to compare iteratively for hamming filtering (Default: 1000)", default=1000)
    parser_clone_umi.add_argument("--pt", type=str, help="Shuffled UMI file path if already made (Default: None)", default=argparse.SUPPRESS)
    parser_clone_umi.add_argument("--dir", type=str, help="Output directory", default='../out')

    parser_clone_umi.set_defaults(func=cl.umi)

    # pcr_sim(): returns dataframe with simulated pcr product 
    parser_clone_pcrsim = subparsers_clone.add_parser("pcr_sim", help="Simulate PCR product from template and primer sequences", description="Simulate PCR product from template and primer sequences", formatter_class=formatter_class)

    parser_clone_pcrsim.add_argument("--df", type=str, help="Input dataframe or file path containing template and primers", required=True)
    parser_clone_pcrsim.add_argument("--template_col", type=str, help="Column name for template sequence", required=True)
    parser_clone_pcrsim.add_argument("--fwd_bind_col", type=str, help="Column name for forward primer binding region", required=True)
    parser_clone_pcrsim.add_argument("--rev_bind_col", type=str, help="Column name for reverse primer binding region", required=True)

    parser_clone_pcrsim.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_clone_pcrsim.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_pcr_sim.csv')

    parser_clone_pcrsim.add_argument("--fwd_ext_col", type=str, help="Column name for forward primer extension region")
    parser_clone_pcrsim.add_argument("--rev_ext_col", type=str, help="Column name for reverse primer extension region")
    parser_clone_pcrsim.add_argument("--product_col", type=str, default="PCR Product", help="Column name for output PCR product")
    parser_clone_pcrsim.add_argument("--no_primer_in_product", dest="primer_in_product", action="store_false",  default=True, help="don't include primer binding and extension regions in PCR product (Default: True = include primers)")
 
    parser_clone_pcrsim.set_defaults(func=cl.pcr_sim)

    # off_targets(): Find off-target sequences for a list of sequences using pairwise alignment.
    parser_clone_off_targets = subparsers_clone.add_parser("off_targets", help="Find off-target sequences for a list of sequences of the same length using pairwise alignment", description="Find off-target sequences for a list of sequences of the same length using pairwise alignment", formatter_class=formatter_class)

    parser_clone_off_targets.add_argument("--df", type=str, help="Input file path containing sequences", required=True)
    parser_clone_off_targets.add_argument("--col", type=str, help="Column name containing sequences to align", required=True)

    parser_clone_off_targets.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_clone_off_targets.add_argument("--ckpt", type=int, help="Checkpoint interval for saving progress (Default: 100)", default=100)
    parser_clone_off_targets.add_argument("--match_score", type=float, help="Match score for pairwise alignment", default=argparse.SUPPRESS)
    parser_clone_off_targets.add_argument("--mismatch_score", type=float, help="Mismatch score for pairwise alignment", default=argparse.SUPPRESS)
    parser_clone_off_targets.add_argument("--open_gap_score", type=float, help="Open gap score for pairwise alignment", default=argparse.SUPPRESS)
    parser_clone_off_targets.add_argument("--extend_gap_score", type=float, help="Extend gap score for pairwise alignment", default=argparse.SUPPRESS)

    parser_clone_off_targets.set_defaults(func=cl.off_targets)

    '''
    edms.bio.transfect
    - PE3(): generates PE3 transfection plan for HEK293T cells (Default: 96-well plate in triplicate using L2000)
    - virus(): generates transfection plan for virus production from HEK293T cells (Default: 6-well plate using L3000)
    '''
    parser_transfect = subparsers.add_parser("transfect", help="Transfection", description="Transfection", formatter_class=formatter_class)
    subparsers_transfect = parser_transfect.add_subparsers()

    # PE3(): generates PE3 transfection plan for HEK293T cells (Default: 96-well plate in triplicate using L2000)
    parser_transfect_PE3 = subparsers_transfect.add_parser("PE3", help="Plan PE3 transfection", description="Plan PE3 transfection", formatter_class=formatter_class)
    
    parser_transfect_PE3.add_argument("--plasmids", type=str, help="Path to plasmids file", required=True)
    parser_transfect_PE3.add_argument("--epegRNAs", type=str, help="Path to epegRNAs file", required=True)
    parser_transfect_PE3.add_argument("--ngRNAs", type=str, help="Path to ngRNAs file", required=True)

    parser_transfect_PE3.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_transfect_PE3.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_transfect_PE3.csv')

    parser_transfect_PE3.add_argument("--pegRNA_number_col", type=str, default="pegRNA_number", help="Column name for pegRNA number")
    parser_transfect_PE3.add_argument("--epegRNAs_name_col", type=str, default="Name", help="Column name for epegRNA name")
    parser_transfect_PE3.add_argument("--ngRNAs_name_col", type=str, default="Name", help="Column name for ngRNA name")
    parser_transfect_PE3.add_argument("--plasmid_col", type=str, default="Plasmid", help="Column name for plasmid name")
    parser_transfect_PE3.add_argument("--description_col", type=str, default="Description", help="Column name for plasmid description")
    parser_transfect_PE3.add_argument("--colony_col", type=str, default="Colony", help="Column name for colony name")
    parser_transfect_PE3.add_argument("--ng_uL_col", type=str, default="ng/uL", help="Column name for ng/uL concentration")
    parser_transfect_PE3.add_argument("--PE_plasmid", type=str, default="pMUZ86.7", help="Name of PE plasmid to search for")
    parser_transfect_PE3.add_argument("--reps", type=int, default=3, help="Number of replicates")
    parser_transfect_PE3.add_argument("--mm_x", type=float, default=1.1, help="Master mix multiplier")
    parser_transfect_PE3.add_argument("--epegRNA_ng", type=int, default=66, help="ng of epegRNA per well")
    parser_transfect_PE3.add_argument("--ngRNA_ng", type=int, default=22, help="ng of ngRNA per well")
    parser_transfect_PE3.add_argument("--PE_ng", type=int, default=200, help="ng of PE plasmid per well")
    parser_transfect_PE3.add_argument("--well_uL", type=int, default=10, help="Total uL per well")

    parser_transfect_PE3.set_defaults(func=tf.PE3)

    # virus(): generates transfection plan for virus production from HEK293T cells (Default: 6-well plate using L3000)
    parser_transfect_virus = subparsers_transfect.add_parser("virus", help="Plan virus transfection", description="Plan virus transfection", formatter_class=formatter_class)

    parser_transfect_virus.add_argument("--plasmids", type=str, help="Path to plasmids file", required=True)

    parser_transfect_virus.add_argument("--dir", type=str, help="Output directory", default='../out')
    parser_transfect_virus.add_argument("--file", type=str, help="Output file name", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_transfect_virus.csv')

    parser_transfect_virus.add_argument("--plasmid_col", type=str, default="Plasmid", help="Column name for plasmid name")
    parser_transfect_virus.add_argument("--description_col", type=str, default="Description", help="Column name for plasmid description")
    parser_transfect_virus.add_argument("--colony_col", type=str, default="Colony", help="Column name for colony")
    parser_transfect_virus.add_argument("--ng_uL_col", type=str, default="ng/uL", help="Column name for ng/uL concentration")
    parser_transfect_virus.add_argument("--VSVG_plasmid", type=str, default="pMUZ26.6", help="Name of VSVG plasmid")
    parser_transfect_virus.add_argument("--GagPol_plasmid", type=str, default="pMUZ26.7", help="Name of GagPol plasmid")
    parser_transfect_virus.add_argument("--reps", type=int, default=1, help="Number of replicates")
    parser_transfect_virus.add_argument("--mm_x", type=float, default=1.1, help="Master mix multiplier")
    parser_transfect_virus.add_argument("--VSVG_ng", type=int, default=750, help="VSVG ng per well")
    parser_transfect_virus.add_argument("--GagPol_ng", type=int, default=1500, help="GagPol ng per well")
    parser_transfect_virus.add_argument("--transfer_ng", type=int, default=750, help="Transfer plasmid ng per well")
    parser_transfect_virus.add_argument("--well_uL", type=int, default=500, help="Total uL per well")

    parser_transfect_virus.set_defaults(func=tf.virus)

    '''
    edms.bio.qPCR:
    - ddCq(): computes ΔΔCq mean and error for all samples holding target pairs constant
    '''
    # ddCq(): computes ΔΔCq mean and error for all samples holding target pairs constant
    parser_ddcq = subparsers.add_parser("ddCq", help="Compute ΔΔCq values for RT-qPCR data", description="Compute ΔΔCq values for RT-qPCR data", formatter_class=formatter_class)
    
    parser_ddcq.add_argument("--data", type=str, help="Input Cq file from CFX instrument",required=True)

    parser_ddcq.add_argument("--dir", type=str, help="Output directory",default='../out')
    parser_ddcq.add_argument("--file", type=str, help="Output file name",default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_qPCR_ddCq.csv')

    parser_ddcq.add_argument("--sample_col", type=str, default="Sample", help="Column name for sample ID")
    parser_ddcq.add_argument("--target_col", type=str, default="Target", help="Column name for target gene ID")
    parser_ddcq.add_argument("--Cq_col", type=str, default="Cq", help="Column name for Cq values")

    parser_ddcq.set_defaults(func=qPCR.ddCq)

    '''
    edms.bio.fastq:
    - savemoney(): create savemoney samples.csv for nanopore mixed WPS
    
    - revcom_fastqs(): write reverse complement of fastqs to a new directory
    - unzip_fastqs(): Unzip gzipped fastqs and write to a new directory
    - comb_fastqs(): Combines one or more (un)compressed fastqs files into a single (un)compressed fastq file
    
    - genotyping(): quantify edit outcomes workflow
    - abundances(): quantify desired edits count & fraction per sample
    
    - count_motif(): returns a dataframe with the sequence motif location with mismatches per read for every fastq file in a directory
    - plot_motif(): generate plots highlighting motif mismatches, locations, and sequences
    - plot_alignments(): generate line & distribution plots from fastq alignments dictionary
    - count_region(): align read region from fastq directory to the annotated library with mismatches; plot and return fastq alignments dictionary
    - count_alignments(): align reads from fastq directory to annotated library with mismatches; plot and return fastq alignments dictionary
    - plot_paired(): generate stacked bar plots from paired_regions() dataframe
    - paired_regions(): quantify, plot, & return (un)paired regions that aligned to the annotated library
    
    - count_signatures(): generate signatures from fastq read region alignments to WT sequence; count signatures, plot and return fastq signatures dataframe
    - editing_per_library(): Determine editing relative library abundance
    
    - extract_umis(): extract UMIs using umi_tools
    - trim_motifs(): trimming motifs with cutadapt
    - make_sams(): generates alignments saved as a SAM files using bowtie2
    - make_bams(): converts SAM files to BAM files using samtools
    - bam_umi_tags(): copy UMI in read ID to RX tag in BAM files using fgbio
    - group_umis(): group BAM files by UMI using fgbio
    - consensus_umis(): generate consensus sequences from grouped UMIs using fgbio
    - bam_to_fastq(): convert BAM files to FASTQ files using samtools

    - cat(): create categorical graphs
    - stack(): create stacked bar plot
    - vol(): create volcano plot
    - torn(): create tornado plot
    - corr(): create correlation plot
    - heat(): create heatmap plot
    '''
    parser_fastq = subparsers.add_parser("fastq", help="FASTQ files", description="FASTQ files", formatter_class=formatter_class)
    subparsers_fastq = parser_fastq.add_subparsers()

    parser_fastq_savemoney = subparsers_fastq.add_parser("savemoney", help="Create savemoney samples.csv for nanopore mixed WPS", description="Create savemoney samples.csv for nanopore mixed WPS", formatter_class=formatter_class)
    
    parser_fastq_revcom = subparsers_fastq.add_parser("revcom", help="Reverse complement all FASTQ files in a directory", description="Reverse complement all FASTQ files in a directory", formatter_class=formatter_class)
    parser_fastq_unzip = subparsers_fastq.add_parser("unzip", help="Unzip gzipped FASTQ files to a new directory", description="Unzip gzipped FASTQ files to a new directory", formatter_class=formatter_class)
    parser_fastq_comb = subparsers_fastq.add_parser("comb", help="Combine multiple FASTQ files into a single FASTQ (.fastq or .fastq.gz)", description="Combine multiple FASTQ files into a single FASTQ (.fastq or .fastq.gz)", formatter_class=formatter_class)
    
    parser_fastq_genotyping = subparsers_fastq.add_parser("genotyping", help="Quantify edit outcomes workflow", description="Quantify edit outcomes workflow", formatter_class=formatter_class)
    parser_fastq_abundances = subparsers_fastq.add_parser("abundances", help="Quantify edit outcomes count & fraction per sample", description="Quantify edit outcomes count & fraction per sample", formatter_class=formatter_class)
    
    parser_fastq_count_motif = subparsers_fastq.add_parser("count_motif", help="Count motif occurrences in FASTQ files", description="Count motif occurrences in FASTQ files", formatter_class=formatter_class)
    parser_fastq_plot_motif = subparsers_fastq.add_parser("plot_motif", help="Plot motif occurrences from FASTQ files", description="Plot motif occurrences from FASTQ files", formatter_class=formatter_class)
    parser_fastq_plot_alignments = subparsers_fastq.add_parser("plot_alignments", help="Plot alignments from FASTQ files", description="Plot alignments from FASTQ files", formatter_class=formatter_class)
    parser_fastq_count_region = subparsers_fastq.add_parser("count_region", help="Count region occurrences in FASTQ files", description="Count region occurrences in FASTQ files", formatter_class=formatter_class)
    parser_fastq_count_alignments = subparsers_fastq.add_parser("count_alignments", help="Count alignments in FASTQ files", description="Count alignments in FASTQ files", formatter_class=formatter_class)
    parser_fastq_plot_paired = subparsers_fastq.add_parser("plot_paired", help="Plot paired regions from FASTQ files", description="Plot paired regions from FASTQ files", formatter_class=formatter_class)
    parser_fastq_paired_regions = subparsers_fastq.add_parser("paired_regions", help="Extract paired regions from FASTQ files", description="Extract paired regions from FASTQ files", formatter_class=formatter_class)
    
    parser_fastq_count_signatures = subparsers_fastq.add_parser("count_signatures", help="Generate signatures from fastq read region alignments to WT sequence", description="Generate signatures from fastq read region alignments to WT sequence", formatter_class=formatter_class)
    parser_fastq_editing_per_library = subparsers_fastq.add_parser("editing_per_library", help="Determine editing relative library abundance", description="Determine editing relative library abundance", formatter_class=formatter_class)
    
    parser_fastq_extract_umis = subparsers_fastq.add_parser("extract_umis", help="Extract UMIs using umi_tools", description="Extract UMIs using umi_tools", formatter_class=formatter_class)
    parser_fastq_trim_motifs = subparsers_fastq.add_parser("trim_motifs", help="Trim motifs with cutadapt", description="Trim motifs with cutadapt", formatter_class=formatter_class)
    parser_fastq_make_sams = subparsers_fastq.add_parser("make_sams", help="Generate alignments saved as SAM files using bowtie2", description="Generate alignments saved as SAM files using bowtie2", formatter_class=formatter_class)
    parser_fastq_make_bams = subparsers_fastq.add_parser("make_bams", help="Convert SAM files to BAM files using samtools", description="Convert SAM files to BAM files using samtools", formatter_class=formatter_class)
    parser_fastq_bam_umi_tags = subparsers_fastq.add_parser("bam_umi_tags", help="Copy UMI in read ID to RX tag in BAM files using fgbio", description="Copy UMI in read ID to RX tag in BAM files using fgbio", formatter_class=formatter_class)
    parser_fastq_group_umis = subparsers_fastq.add_parser("group_umis", help="Group BAM files by UMI using fgbio", description="Group BAM files by UMI using fgbio", formatter_class=formatter_class)
    parser_fastq_consensus_umis = subparsers_fastq.add_parser("consensus_umis", help="Generate consensus sequences from grouped UMIs using fgbio", description="Generate consensus sequences from grouped UMIs using fgbio", formatter_class=formatter_class)
    parser_fastq_bam_to_fastq = subparsers_fastq.add_parser("bam_to_fastq", help="Convert BAM files to FASTQ files using samtools", description="Convert BAM files to FASTQ files using samtools", formatter_class=formatter_class)
    
    parser_fastq_cat = subparsers_fastq.add_parser("cat", help="Create categorical graphs", description="Create categorical graphs", formatter_class=formatter_class)
    parser_fastq_stack = subparsers_fastq.add_parser("stack", help="Create stacked bar plot", description="Create stacked bar plot", formatter_class=formatter_class)
    parser_fastq_vol = subparsers_fastq.add_parser("vol", help="Create volcano plot", description="Create volcano plot", formatter_class=formatter_class)
    parser_fastq_torn = subparsers_fastq.add_parser("torn", help="Create tornado plot", description="Create tornado plot", formatter_class=formatter_class)
    parser_fastq_corr = subparsers_fastq.add_parser("corr", help="Create correlation plot", description="Create correlation plot", formatter_class=formatter_class)
    parser_fastq_heat = subparsers_fastq.add_parser("heat", help="Create heatmap plot", description="Create heatmap plot", formatter_class=formatter_class)

    # savemoney():
    parser_fastq_savemoney.add_argument("--fastq_dir", type=str, help="Path to fastq directory (contains .fastq files, not .fastq.gz files; Default: './fastq')", default='./fastq')
    parser_fastq_savemoney.add_argument("--fasta_dir", type=str, help="Path to fasta directory (contains .fasta files; Default: './fasta')", default='./fasta')
    parser_fastq_savemoney.add_argument("--pt", type=str, help="Working directory when running savemoney (Default: '.' = current directory)", default=".")
    parser_fastq_savemoney.add_argument("--out_dir", type=str, help="Path to output directory (Default: '.' = current directory)", default='.')
    parser_fastq_savemoney.add_argument("--out_file", type=str, help="Name of output file (Default: 'samples.csv')", default='samples.csv')

    # Add common arguments: revcom_fastqs(), unzip_fastqs(), comb_fastqs(), and genotyping()
    for parser_fastq_common in [parser_fastq_revcom,parser_fastq_unzip,parser_fastq_comb,parser_fastq_genotyping]:
        parser_fastq_common.add_argument("--in_dir", type=str, help="Input directory containing FASTQ files",default='.')
        parser_fastq_common.add_argument("--out_dir", type=str, help="Output directory",default = f'../out')

    # Add specific arguments: comb_fastqs()
    parser_fastq_comb.add_argument("--out_file", type=str, help="Name of output FASTQ file (.fastq or .fastq.gz)", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_comb.fastq.gz')

    # Add specific arguments: genotyping()
    parser_fastq_genotyping.add_argument("--out_file_prefix", type=str, help="Name of output file prefix", default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    
    parser_fastq_genotyping.add_argument("--config_key", type=str, help="Config file key (FWD primer-REV primer) for sequence [flank5(genotype region)flank3 or genotype region only] & res (first amino acid in genotype region)", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--sequence", type=str, help="Formatted sequence: flank5(genotype region)flank3 or genotype region only", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--res", type=int, help="First amino acid number in genotype region", default=argparse.SUPPRESS)
    
    parser_fastq_genotyping.add_argument("--qall", type=int, help="Minimum Phred quality score for all bases", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--qtrim", type=int, help="Phred quality threshold for end trimming", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--qavg", type=int, help="Minimum average Phred quality score", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--qmask", type=int, help="Phred quality threshold for masking to N", default=argparse.SUPPRESS)

    parser_fastq_genotyping.add_argument("--save", action="store_true", help="Save read statistics and genotypes files", dest="save", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--no_save", action="store_false", help="Don't save read statistics and genotypes files", dest="save", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--masks", action="store_true", help="Include masked sequence and translation",default=False)
    parser_fastq_genotyping.add_argument("--keepX", action="store_true", help="Keep unknown translation (X) in output", default=False)

    parser_fastq_genotyping.add_argument("--match_score", type=float, help="Match score for pairwise alignment", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--mismatch_score", type=float, help="Mismatch score for pairwise alignment", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--open_gap_score", type=float, help="Open gap score for pairwise alignment", default=argparse.SUPPRESS)
    parser_fastq_genotyping.add_argument("--extend_gap_score", type=float, help="Extend gap score for pairwise alignment", default=argparse.SUPPRESS)
    
    # abundances():
    parser_fastq_abundances.add_argument("--df", help="Input file with sample, edit, count, & fraction information", required=True)
    parser_fastq_abundances.add_argument("--desired_edits", nargs="+", help="List of desired edits to isolate (space-separated)",required=True)

    parser_fastq_abundances.add_argument("--sample_col", default="sample", help="Column for sample ID (Default: 'sample')")
    parser_fastq_abundances.add_argument("--edit_col", default="edit", help="Column for edit identifier (Default: 'edit')")
    parser_fastq_abundances.add_argument("--count_col", default="count", help="Column for edit count (Default: 'count')")
    parser_fastq_abundances.add_argument("--fraction_col", default="fraction", help="Column for edit fraction (Default: 'fraction')")
    parser_fastq_abundances.add_argument("--combinations", default=1, help="Maximum # of desired edit combinations to search for (Default: 1 => single edits)")

    # count_motif():
    parser_fastq_count_motif.add_argument("--fastq_dir", help="Path to directory containing FASTQ files", required=True)
    parser_fastq_count_motif.add_argument("--pattern", help="Motif sequence pattern to search for", required=True)
    parser_fastq_count_motif.add_argument("--out_dir", help="Output directory to save results", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')

    parser_fastq_count_motif.add_argument("--motif", default="motif", help="Name of the motif (Default: 'motif')")
    parser_fastq_count_motif.add_argument("--max_distance", type=int, default=0, help="Maximum Levenshtein distance allowed (i.e., # of mismatches, Default: 0)")
    parser_fastq_count_motif.add_argument("--max_reads", type=int, default=0, help="Maximum # of reads to process per file")
    parser_fastq_count_motif.add_argument("--meta", type=str, help="Optional path to metadata CSV/TSV file with 'fastq_file' column")

    # plot_motif():
    parser_fastq_plot_motif.add_argument("--df", help="Path to count_motif() output file", required=True)

    parser_fastq_plot_motif.add_argument("--out_dir", type=str, help="Directory to save plots", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_plot_motif.add_argument("--plot_suf", type=str, default=".pdf", help="Plot file suffix (e.g. .pdf, .png)")
    parser_fastq_plot_motif.add_argument("--numeric", choices=["count", "fraction"], default="count", help="Numeric column to use for plotting (Default: 'count')")
    parser_fastq_plot_motif.add_argument("--id_col", default="fastq_file", help="Column used for sample ID (Default: 'fastq_file')")
    parser_fastq_plot_motif.add_argument("--id_axis", default="fastq", help="Label to use on the plot axis (Default: 'fastq')")
    parser_fastq_plot_motif.add_argument("--stack_figsize", type=parse_tuple_int, default=(7, 3), help="Stacked plot figure size formatted as 'width,height'")
    parser_fastq_plot_motif.add_argument("--heat_figsize", type=parse_tuple_int, help="Heatmap figure size formateed as 'width,height'")
    parser_fastq_plot_motif.add_argument("--cutoff_frac", type=float, default=0.01, help="Minimum fraction to include in y-axis (Default: 0.01)")

    # plot_alignments():
    parser_fastq_plot_alignments.add_argument("--fastq_alignments", help="Directory with the fastq alignments dictionary", required=True)
    parser_fastq_plot_alignments.add_argument("--align_col", help="Align column name in the annotated library reference file", required=True)
    parser_fastq_plot_alignments.add_argument("--id_col", help="ID column name in the annotated library reference file", required=True)

    parser_fastq_plot_alignments.add_argument("--out_dir", help="Output directory for plots", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_plot_alignments.add_argument("--plot_suf", default=".pdf", help="Plot file suffix (Default: .pdf)")
    parser_fastq_plot_alignments.add_argument("--show", action="store_true", help="Display plots interactively",default=False)
    
    # count_region():
    parser_fastq_count_region.add_argument("--df_ref", help="Annotated reference library file path", required=True)
    parser_fastq_count_region.add_argument("--align_col", help="Align column name in the annotated reference library", required=True)
    parser_fastq_count_region.add_argument("--id_col", help="ID column name in the annotated reference library", required=True)
    parser_fastq_count_region.add_argument("--fastq_dir", help="Directory containing FASTQ files", required=True)
    parser_fastq_count_region.add_argument("--df_motif5", help="5' motif file path", required=True)
    parser_fastq_count_region.add_argument("--df_motif3", help="3' motif file path", required=True)

    parser_fastq_count_region.add_argument("--out_dir", help="Output directory", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_count_region.add_argument("--fastq_col", help="Fastq column name in the annotated reference library (Default: None)", default=None)
    parser_fastq_count_region.add_argument("--match_score", type=float, default=2, help="Score for matches (Default: 2)")
    parser_fastq_count_region.add_argument("--mismatch_score", type=float, default=-1, help="Score for mismatches (Default: -1)")
    parser_fastq_count_region.add_argument("--open_gap_score", type=float, default=-10, help="Gap opening score (Default: -10)")
    parser_fastq_count_region.add_argument("--extend_gap_score", type=float, default=-0.1, help="Gap extension score (Default: -0.1)")
    parser_fastq_count_region.add_argument("--align_dims", type=parse_tuple_int, default=(0, 0), help="Alignment range formatted as 'start,end' (Default: 0,0 = all reads)")
    parser_fastq_count_region.add_argument("--align_ckpt", type=int, default=10000, help="Checkpoint frequency (Default: 10000)")
    parser_fastq_count_region.add_argument("--plot_suf", type=str, help="Plot suffix type (e.g. '.pdf')")
    parser_fastq_count_region.add_argument("--show", action="store_true", help="Display plots interactively", default=False)
    parser_fastq_count_region.add_argument("--exact", action="store_true", help="Perform exact matching only", default=False)
    
    # count_alignments():
    parser_fastq_count_alignments.add_argument("--df_ref", help="Annotated reference library file path", required=True)
    parser_fastq_count_alignments.add_argument("--align_col", help="Align column name in the annotated reference library", required=True)
    parser_fastq_count_alignments.add_argument("--id_col", help="ID column name in the annotated reference library", required=True)
    parser_fastq_count_alignments.add_argument("--fastq_dir", help="Directory containing FASTQ files", required=True)

    parser_fastq_count_alignments.add_argument("--out_dir", help="Output directory", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_count_alignments.add_argument("--fastq_col", help="Fastq column name in the annotated reference library (Default: None)", default=None)
    parser_fastq_count_alignments.add_argument("--match_score", type=float, default=2, help="Match score (Default: 2)")
    parser_fastq_count_alignments.add_argument("--mismatch_score", type=float, default=-1, help="Mismatch penalty (Default: -1)")
    parser_fastq_count_alignments.add_argument("--open_gap_score", type=float, default=-10, help="Gap open penalty (Default: -10)")
    parser_fastq_count_alignments.add_argument("--extend_gap_score", type=float, default=-0.1, help="Gap extension penalty (Default: -0.1)")
    parser_fastq_count_alignments.add_argument("--align_dims", type=parse_tuple_int, default=(0, 0), help="Alignment range as 'start,end' (Default: 0,0 = all reads)")
    parser_fastq_count_alignments.add_argument("--align_ckpt", type=int, default=10000, help="Checkpoint frequency for saving alignment progress")
    parser_fastq_count_alignments.add_argument("--plot_suf", type=str, help="Plot file suffix (e.g. .pdf, .png)")
    parser_fastq_count_alignments.add_argument("--show", action="store_true", help="Show plots interactively")
    parser_fastq_count_alignments.add_argument("--exact", action="store_true", help="Perform exact matching only", default=False)
    
    # plot_paired():
    parser_fastq_plot_paired.add_argument("--df", help="Paired region file path", required=True)
    parser_fastq_plot_paired.add_argument("--title", help="Plot title and output filename (without extension)", required=True)
    
    parser_fastq_plot_paired.add_argument("--out_dir", help="Output directory", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_plot_paired.add_argument("--id_col", default="ID", help="Column name for ID (Default: 'ID')")
    parser_fastq_plot_paired.add_argument("--desired_col", default="desired", help="Column name for desired sequences (Default: 'desired')")
    parser_fastq_plot_paired.add_argument("--y", default="count", help="y axis for plots (Default: 'count'; Options: 'count' & 'fraction')")
    parser_fastq_plot_paired.add_argument("--plot_suf", default=".pdf", help="Plot file suffix (e.g. .pdf or .png)")
    parser_fastq_plot_paired.add_argument("--show", action="store_true", help="Display plots interactively")

    # paired_regions():
    parser_fastq_paired_regions.add_argument("--meta_dir", help="Directory containing meta files", required=True)
    parser_fastq_paired_regions.add_argument("--region1_dir", help="Directory with region 1 alignment files", required=True)
    parser_fastq_paired_regions.add_argument("--region2_dir", help="Directory with region 2 alignment files", required=True)

    parser_fastq_paired_regions.add_argument("--out_dir", help="Output directory", default=f'../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}')
    parser_fastq_paired_regions.add_argument("--id_col", default="ID", help="Column name for unique identifiers (Default: 'ID')")
    parser_fastq_paired_regions.add_argument("--desired_col", default="desired", help="Column name for desired sequences (Default: 'desired')")
    parser_fastq_paired_regions.add_argument("--region1_alignment_col", default="r1_alignment", help="Column name for region 1 alignment data")
    parser_fastq_paired_regions.add_argument("--region2_alignment_col", default="r2_alignment", help="Column name for region 2 alignment data")
    parser_fastq_paired_regions.add_argument("--reads_aligned_col", default="reads_aligned", help="Column name for aligned reads (Default: 'reads_aligned')")
    parser_fastq_paired_regions.add_argument("--reads_processed_col", default="reads_processed", help="Column name for processed reads (Default: 'reads_processed')")
    parser_fastq_paired_regions.add_argument("--plot_suf", default=".pdf", help="Plot file suffix (e.g., .pdf, .png)")
    parser_fastq_paired_regions.add_argument("--show", action="store_true", help="Display plots interactively")
    parser_fastq_paired_regions.add_argument("--return_dc", action="store_true", help="Return paired/unpaired dataframe")

    # count_signatures():
    parser_fastq_count_signatures.add_argument("--df_ref", help="Annotated reference library file path", required=True)
    parser_fastq_count_signatures.add_argument("--fastq_dir", help="Directory containing FASTQ files", required=True)
    
    parser_fastq_count_signatures_group = parser_fastq_count_signatures.add_mutually_exclusive_group(required=True)
    parser_fastq_count_signatures_group.add_argument("--config_key", type=str, help="[Required (Option 1)] Config file key (FWD primer_REV primer) with 'motif5' and 'motif3'", default=argparse.SUPPRESS)
    parser_fastq_count_signatures_group.add_argument("--in_file", help="[Required (Option 2)] Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence,index (column names required)", default=argparse.SUPPRESS)
    parser_fastq_count_signatures_group.add_argument("--sequence", help="[Required (Option 3)] Target sequence; retrieved from config_key or in_file if not provided")
    
    parser_fastq_count_signatures.add_argument("--n_extra_nt", type=int, help="Number of extra nucleotide differences allowed for Signature match (Default: 0)", default=0)
    parser_fastq_count_signatures.add_argument("--df_motif5", help="5' motif file path", default=argparse.SUPPRESS)
    parser_fastq_count_signatures.add_argument("--df_motif3", help="3' motif file path", default=argparse.SUPPRESS)
    parser_fastq_count_signatures.add_argument("--meta", help="Meta file path", default=argparse.SUPPRESS)
    parser_fastq_count_signatures.add_argument("--out_dir", help="Output directory", default='../out/')
    parser_fastq_count_signatures.add_argument("--out_file", help="Output filename", default=f"{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_count_signatures.csv")
    parser_fastq_count_signatures.add_argument("--signature_col", help="Signature column name in the annotated reference library (Default: 'Signature')", default='Signature')
    parser_fastq_count_signatures.add_argument("--id_col", help="ID column name in the annotated reference library (Default: 'ID')", default='ID')
    parser_fastq_count_signatures.add_argument("--edit_col", help="Edit column name in the annotated reference library (Default: 'Edit')", default='Edit')
    parser_fastq_count_signatures.add_argument("--fastq_col", help="Fastq column name in the annotated reference library (Default: None)", default=None)
    parser_fastq_count_signatures.add_argument("--match_score", type=float, default=2, help="Score for matches (Default: 2)")
    parser_fastq_count_signatures.add_argument("--mismatch_score", type=float, default=-1, help="Score for mismatches (Default: -1)")
    parser_fastq_count_signatures.add_argument("--open_gap_score", type=float, default=-10, help="Gap opening score (Default: -10)")
    parser_fastq_count_signatures.add_argument("--extend_gap_score", type=float, default=-0.1, help="Gap extension score (Default: -0.1)")
    parser_fastq_count_signatures.add_argument("--align_dims", type=parse_tuple_int, default=(0, 0), help="Alignment range formatted as 'start,end' (Default: 0,0 = all reads)")
    parser_fastq_count_signatures.add_argument("--align_ckpt", type=int, default=10000, help="Checkpoint frequency (Default: 10000)")
    parser_fastq_count_signatures.add_argument("--save_alignments", action="store_true", help="Save alignments (Default: False, save memory)", default=False)
    parser_fastq_count_signatures.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)
    parser_fastq_count_signatures.add_argument("--plot_suf", type=str, help="Plot suffix type (Default: 'pdf')", default='.pdf')
    parser_fastq_count_signatures.add_argument("--show", action="store_true", help="Display plots interactively", default=False)
    
    # editing_per_library():
    parser_fastq_editing_per_library.add_argument("--edit_dc", help="Path to directory with edit outcomes files", required=True)
    parser_fastq_editing_per_library.add_argument("--paired_regions_dc", help="Path to directory with paired regions files", required=True)
    parser_fastq_editing_per_library.add_argument("--fastq_ids", help="Path to file containing fastq IDs for 'genotyping' and 'paired_regions'", required=True)

    parser_fastq_editing_per_library.add_argument("--out_dir", type=str, help="Output directory to save results (Default: ../out/date_time)", default=f"../out/{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}")
    parser_fastq_editing_per_library.add_argument("--count", default="count", help="Column to use for epeg-ngRNA counts (Default: 'count')")
    parser_fastq_editing_per_library.add_argument("--psuedocount", type=int, default=1, help="Pseudocount to add to all counts (Default: 1)")

    # extract_umis():
    parser_fastq_extract_umis.add_argument("--fastq_dir", help="Directory containing FASTQ files", required=True)

    parser_fastq_extract_umis.add_argument("--out_dir", help="Output directory (Default: ./extract_umis)", default=f'./extract_umis')
    parser_fastq_extract_umis.add_argument("--bc_pattern", help="UMI barcode pattern (Default: NNNNNNNNNNNNNNNN)", default="NNNNNNNNNNNNNNNN")
    parser_fastq_extract_umis.add_argument("--env", help="Conda environment with umi_tools installed (Default: umi_tools)", default="umi_tools")

    # trim_motifs():
    parser_fastq_trim_motifs.add_argument("--fastq_dir", help="Directory containing FASTQ files (with UMIs extracted)", required=True)

    parser_fastq_trim_motifs.add_argument("--out_dir", help="Output directory (Default: ./trim_motifs)", default=f'./trim_motifs')

    parser_fastq_trim_motifs.add_argument("--config_key", type=str, help="[Required (Option 1)] Config file key (FWD primer_REV primer) with 'motif5' and 'motif3'", default=argparse.SUPPRESS)
    parser_fastq_trim_motifs.add_argument("--in_file", help="[Required (Option 2)] Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)", default=argparse.SUPPRESS)
    parser_fastq_trim_motifs.add_argument("--motif5", help="[Required (Option 3 and/or)] 5' motif sequence to trim", default=argparse.SUPPRESS)
    parser_fastq_trim_motifs.add_argument("--motif3", help="[Required (Option 3 and/or)] 3' motif sequence to trim", default=argparse.SUPPRESS)

    parser_fastq_trim_motifs.add_argument("--motif_length", type=int, help="Trim 'in_file' motifs to this length (Default: 21)", default=21)
    parser_fastq_trim_motifs.add_argument("--error_rate", type=float, help="Maximum error rate allowed in each motif (Default: 0.1 = 10%%)", default=0.1)
    parser_fastq_trim_motifs.add_argument("--env", help="Conda environment with cutadapt installed (Default: umi_tools)", default="umi_tools")

    # make_sams():
    parser_fastq_make_sams.add_argument("--fastq_dir", help="Directory containing FASTQ files", required=True)

    parser_fastq_make_sams.add_argument("--out_dir", help="Output directory (Default: ./make_sams)", default=f'./make_sams')

    parser_fastq_make_sams_group = parser_fastq_make_sams.add_mutually_exclusive_group(required=True)
    parser_fastq_make_sams_group.add_argument("--in_file", help="Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence,index (column names required)", default=argparse.SUPPRESS)
    parser_fastq_make_sams_group.add_argument("--fasta", help="Reference FASTA file for alignment", default=argparse.SUPPRESS)

    parser_fastq_make_sams.add_argument("--sensitivity", choices=["very-sensitive", "sensitive", "fast", "very-fast", "very-sensitive-local", "sensitive-local", "fast-local", "very-fast-local"], 
                                        default="very-sensitive", help="Bowtie2 sensitivity setting (Default: very-sensitive)")

    parser_fastq_make_sams.add_argument("--env", help="Conda environment with bowtie2 installed (Default: umi_tools)", default="umi_tools")

    # make_bams():
    parser_fastq_make_bams.add_argument("--sam_dir", help="Directory containing SAM files", required=True)

    parser_fastq_make_bams.add_argument("--out_dir", help="Output directory (Default: ./make_bams)", default=f'./make_bams')
    parser_fastq_make_bams.add_argument("--env", help="Conda environment with samtools installed (Default: umi_tools)", default="umi_tools")

    # bam_umi_tags():
    parser_fastq_bam_umi_tags.add_argument("--bam_dir", help="Directory containing bam files", required=True)

    parser_fastq_bam_umi_tags.add_argument("--out_dir", help="Output directory (Default: ./bam_umi_tags)", default=f'./bam_umi_tags')
    parser_fastq_bam_umi_tags.add_argument("--env", help="Conda environment with fgbio installed (Default: umi_tools)", default="umi_tools")

    # group_umis():
    parser_fastq_group_umis.add_argument("--bam_dir", help="Directory containing BAM files", required=True)

    parser_fastq_group_umis.add_argument("--out_dir", help="Output directory (Default: ./group_umis)", default=f'./group_umis')
    parser_fastq_group_umis.add_argument("--strategy", choices=["Identical","Edit","Adjacency", "Paired"], help="umi grouping strategy (Default: Adjacency)", default="Adjacency")
    parser_fastq_group_umis.add_argument("--edits", type=int, help="Maximum edit distance to group UMIs (Default: 1)", default=1)
    parser_fastq_group_umis.add_argument("--env", help="Conda environment with fgbio installed (Default: umi_tools)", default="umi_tools")

    # consensus_umis():
    parser_fastq_consensus_umis.add_argument("--bam_dir", help="Directory containing grouped BAM files", required=True)

    parser_fastq_consensus_umis.add_argument("--out_dir", help="Output directory (Default: ./consensus_umis)", default=f'./consensus_umis')
    parser_fastq_consensus_umis.add_argument("--min_reads", type=int, help="Minimum reads per UMI to call consensus (Default: 1)", default=1)
    parser_fastq_consensus_umis.add_argument("--env", help="Conda environment with fgbio installed (Default: umi_tools)", default="umi_tools")

    # bam_to_fastq():
    parser_fastq_bam_to_fastq.add_argument("--bam_dir", help="Directory containing BAM files", required=True)
    parser_fastq_bam_to_fastq.add_argument("--out_dir", help="Output directory (Default: ./bam_to_fastq)", default=f'./bam_to_fastq')
    parser_fastq_bam_to_fastq.add_argument("--env", help="Conda environment with samtools installed (Default: umi_tools)", default="umi_tools")
    
    # cat():
    add_common_plot_cat_args(parser_fastq_cat, fastq_parser=True)

    # stack():
    add_common_plot_stack_args(parser_fastq_stack, fastq_parser=True)

    # vol():
    add_common_plot_vol_args(parser_fastq_vol, fastq_parser=True)
    
    # torn():
    add_common_plot_scat_args(parser_fastq_torn, fastq_torn_parser=True)

    # corr():
    add_common_plot_scat_args(parser_fastq_corr, fastq_corr_parser=True)

    # heat():
    add_common_plot_heat_args(parser_fastq_heat, fastq_parser=True)
    
    # Set defaults
    parser_fastq_savemoney.set_defaults(func=fq.savemoney)
    parser_fastq_revcom.set_defaults(func=fq.revcom_fastqs)
    parser_fastq_unzip.set_defaults(func=fq.unzip_fastqs)
    parser_fastq_comb.set_defaults(func=fq.comb_fastqs)
    parser_fastq_genotyping.set_defaults(func=fq.genotyping)
    parser_fastq_abundances.set_defaults(func=fq.abundances)
    parser_fastq_count_motif.set_defaults(func=fq.count_motif)
    parser_fastq_plot_motif.set_defaults(func=fq.plot_motif)
    parser_fastq_plot_alignments.set_defaults(func=fq.plot_alignments)
    parser_fastq_count_region.set_defaults(func=fq.count_region)
    parser_fastq_count_alignments.set_defaults(func=fq.count_alignments)
    parser_fastq_plot_paired.set_defaults(func=fq.plot_paired)
    parser_fastq_paired_regions.set_defaults(func=fq.paired_regions)
    parser_fastq_count_signatures.set_defaults(func=fq.count_signatures)
    parser_fastq_editing_per_library.set_defaults(func=fq.editing_per_library)
    parser_fastq_extract_umis.set_defaults(func=fq.extract_umis)
    parser_fastq_trim_motifs.set_defaults(func=fq.trim_motifs)
    parser_fastq_make_sams.set_defaults(func=fq.make_sams)
    parser_fastq_make_bams.set_defaults(func=fq.make_bams)
    parser_fastq_bam_umi_tags.set_defaults(func=fq.bam_umi_tags)
    parser_fastq_group_umis.set_defaults(func=fq.group_umis)
    parser_fastq_consensus_umis.set_defaults(func=fq.consensus_umis)
    parser_fastq_bam_to_fastq.set_defaults(func=fq.bam_to_fastq)
    parser_fastq_cat.set_defaults(func=fq.cat)
    parser_fastq_stack.set_defaults(func=fq.stack)
    parser_fastq_vol.set_defaults(func=fq.vol)
    parser_fastq_torn.set_defaults(func=fq.torn)
    parser_fastq_corr.set_defaults(func=fq.corr)

    '''
    edms.bio.pe:
    - prime_designer(): Execute PrimeDesign saturation mutagenesis (EDMS version)
    - pilot_screen(): Create pilot screen for EDMS
    - epegRNA_linkers(): Generate epegRNA linkers between PBS and 3' hairpin motif & finish annotations
    - merge(): rejoins epeg/ngRNAs & creates ngRNA_groups
    - sensor_designer(): design pegRNA sensors
    - pegRNA_outcome(): confirm that pegRNAs should create the predicted edits
    - pegRNA_signature(): create signatures for pegRNA outcomes using alignments
    - epegRNA_fastas(): generate FASTA files representing epegRNAs cloned into a linearized vector
    '''
    parser_pe = subparsers.add_parser("pe", help="Prime Editing", formatter_class=formatter_class)
    subparsers_pe = parser_pe.add_subparsers()

    parser_pe_prime_designer = subparsers_pe.add_parser("prime_designer", help="Execute PrimeDesign saturation mutagenesis (EDMS version)", description="Execute PrimeDesign saturation mutagenesis (EDMS version)", formatter_class=formatter_class)
    parser_pe_pilot_screen = subparsers_pe.add_parser("pilot_screen", help="Determine pilot screen for EDMS", description="Determine pilot screen for EDMS", formatter_class=formatter_class)
    parser_pe_epegRNA_linkers = subparsers_pe.add_parser("epegRNA_linkers", help="Generate epegRNA linkers between PBS and 3' hairpin motif", description="Generate epegRNA linkers between PBS and 3' hairpin motif", formatter_class=formatter_class)
    parser_pe_merge = subparsers_pe.add_parser("merge", help="rejoins epeg/ngRNAs & creates ngRNA groups", description="rejoins epeg/ngRNAs & creates ngRNA groups", formatter_class=formatter_class)
    parser_pe_sensor_designer = subparsers_pe.add_parser("sensor_designer", help='Design pegRNA sensors', description='Design pegRNA sensors', formatter_class=formatter_class)
    parser_pe_pegRNA_outcome = subparsers_pe.add_parser("pegRNA_outcome", help="Confirm that pegRNAs should create the predicted edit", description="Confirm that pegRNAs should create the predicted edit", formatter_class=formatter_class)
    parser_pe_pegRNA_signature = subparsers_pe.add_parser("pegRNA_signature", help="Create signatures for pegRNA outcomes using alignments", description="Create signatures for pegRNA outcomes using alignments", formatter_class=formatter_class)
    parser_pe_epegRNA_fasta = subparsers_pe.add_parser("epegRNA_fasta", help="Generate FASTA files representing epegRNAs cloned into a linearized vector", description="Generate FASTA files representing epegRNAs cloned into a linearized vector", formatter_class=formatter_class)

    # prime_designer():
    parser_pe_prime_designer.add_argument("--file", type=str, dest='in_file', help="[Required (Option 1)] Input file (.csv or .txt) with sequences for PrimeDesign. Format: target_name,target_sequence,index (Required). See examples below...")
    parser_pe_prime_designer.add_argument("--name", type=str, dest='target_name',help="[Required (Option 2)] Name of the target")
    parser_pe_prime_designer.add_argument("--flank5", type=str, dest='flank5_sequence', help="[Required (Option 2)] 5' flank sequence (in-frame, length divisible by 3)")
    parser_pe_prime_designer.add_argument("--target", type=str, dest='target_sequence', help="[Required (Option 2)] Target sequence (in-frame, length divisible by 3)")
    parser_pe_prime_designer.add_argument("--flank3", type=str, dest='flank3_sequence', help="[Required (Option 2)] 3' flank sequence (in-frame, length divisible by 3)")
    parser_pe_prime_designer.add_argument("--index", type=int, default=1,
                        help="[Required (Option 2)] Index of 1st amino acid or base in target sequence (Default: 1)")
    
    parser_pe_prime_designer.add_argument("--pe_format", type=str, default="NNNNNNNNNNNNNNNNN/NNN[NGG]",
                        help="Prime editing formatting including the spacer, cut index -> /, and protospacer adjacent motif (PAM) -> [PAM] (Default: NNNNNNNNNNNNNNNNN/NNN[NGG]). Warning: Changing pe_format prevents silent mutations from being applied.")
    parser_pe_prime_designer.add_argument("--pbs_length_list", type=int, default=argparse.SUPPRESS, nargs="+",
                        help="List of primer binding site (PBS) lengths for the pegRNA extension")
    parser_pe_prime_designer.add_argument("--rtt_length_list", type=int, default=argparse.SUPPRESS, nargs="+",
                        help="List of reverse transcription template (RTT) lengths for the pegRNA extension")
    parser_pe_prime_designer.add_argument("--nicking_distance_minimum", type=int, default=0,
                        help="Minimum nicking distance for pegRNA designs (Default: 0 nt)")
    parser_pe_prime_designer.add_argument("--nicking_distance_maximum", type=int, default=100,
                        help="Maximum nicking distance for pegRNA designs (Default: 100 nt)")
    parser_pe_prime_designer.add_argument("--filter_c1_extension", action="store_true",
                        help="Filter against pegRNA extensions that start with a C base (Default: False)", default=False)
    parser_pe_prime_designer.add_argument("--no_silent_mutation", dest="silent_mutation", action="store_false",
                        help="Disable silent mutation", default=True)
    parser_pe_prime_designer.add_argument("--genome_wide_design", action="store_true",
                        help="Whether this is a genome-wide pooled design. This option designs a set of pegRNAs per input without ranging PBS and RTT parameters", default=False)
    parser_pe_prime_designer.add_argument("--saturation_mutagenesis", type=str, choices=['aa', 'aa_subs', 'aa_ins', 'aa_dels', 'base'], help="Saturation mutagenesis design with prime editing. The 'aa' option makes all amino acid substitutions ('aa_subs'), +1 amino acid insertions ('aa_ins'), and -1 amino acid deletions ('aa_dels'). The 'base' option makes DNA base changes.", default=None)
    parser_pe_prime_designer.add_argument("--number_of_pegrnas", type=int, default=3,
                        help="Max number of pegRNAs to design (Default: 3)")
    parser_pe_prime_designer.add_argument("--number_of_ngrnas", type=int, default=3,
                        help="Max number of ngRNAs to design (Default: 3)")
    parser_pe_prime_designer.add_argument("--nicking_distance_pooled", type=int, default=75,
                        help="The nicking distance between pegRNAs and ngRNAs for pooled designs. PE3b annotation is priority (PE3b seed -> PE3b non-seed), followed by nicking distance closest to this parameter. (Default: 75 bp)")
    parser_pe_prime_designer.add_argument("--homology_downstream", type=int, default=10,
                        help="Minimum RT extension length downstream of an edit for pegRNA designs (Default: 10 nt)")
    parser_pe_prime_designer.add_argument("--pbs_length_pooled_list", type=int, dest='pbs_length_pooled_list', nargs="+", default=[11,13,15],
                        help="List of PBS lengths to design pegRNAs for pooled design applications (Default: 11,13,15)")
    parser_pe_prime_designer.add_argument("--rtt_max_length_pooled", type=int, default=50,
                        help="Maximum RTT length to design pegRNAs for pooled design applications (Default: 50 nt)")
    parser_pe_prime_designer.add_argument("--scaffold_sequence", type=str, default="GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
                        help="sgRNA scaffold sequence (Default: SpCas9 flip + extend")
    parser_pe_prime_designer.add_argument("--enzymes", type=str, nargs="+", help="list of type IIS RE enzymes (i.e., Esp3I, BsaI, BspMI) to check for in pegRNAs and ngRNAs (Default: ['Esp3I'])", default=['Esp3I'])
    parser_pe_prime_designer.add_argument("--dont_replace", action='store_false',dest='replace', help="Do not replace pegRNAs and remove ngRNAs with RE enzyme sites", default=True)

    # Help message for edms pe prime_designer because input file format can't be captured as block text by Myformatter(RichHelpFormatter):
    if any(["edms" in argv for argv in sys.argv]) and "pe" in sys.argv and "prime_designer" in sys.argv and ("--help" in sys.argv or "-h" in sys.argv):
        parser_pe_prime_designer.print_help()
        rprint("""[red]
Examples:[/red]
  [cyan]--file[/cyan] [dark_magenta]IN_FILE[/dark_magenta]
  [blue]Input file (.csv or .txt) with sequences for PrimeDesign.
  Format: target_name,target_sequence,index

  *** Example saturation_mutagenesis.CSV file *** ---------------------------------------
  |											|
  |	target,ATGTGC(TGTGATGGTATGCCGGCGTAGTAA)TCGTAG,1		                        |
  |											|
  ---------------------------------------------------------------------------------------

  *** Example saturation_mutagenesis.TXT file *** ---------------------------------------
  |											|
  |	target	ATGTGC(TGTGATGGTATGCCGGCGTAGTAA)TCGTAG   1                              |
  |											|
  ---------------------------------------------------------------------------------------

  *** Example not_saturation_mutagenesis.CSV file *** -----------------------------------
  |											|
  |	target_01_substitution,ATGTGCTGTGATGGTAT(G/A)CCGGCGTAGTAATCGTAGC,1		|
  |	target_01_insertion,ATGTGCTGTGATGGTATG(+ATCTCGATGA)CCGGCGTAGTAATCGTAGC,1	|
  |	target_01_deletion,ATGTGCTGTGATGG(-TATGCCG)GCGTAGTAATCGTAGC,1			|
  |											|
  ---------------------------------------------------------------------------------------

  *** Example not_saturation_mutagenesis.TXT file *** -----------------------------------
  |											|
  |	target_01_substitution	ATGTGCTGTGATGGTAT(G/A)CCGGCGTAGTAATCGTAGC   1           |
  |	target_01_insertion	ATGTGCTGTGATGGTATG(+ATCTCGATGA)CCGGCGTAGTAATCGTAGC  1   |
  |	target_01_deletion	ATGTGCTGTGATGG(-TATGCCG)GCGTAGTAATCGTAGC    1           |
  |											|
  ---------------------------------------------------------------------------------------

  *** Formatting different DNA edits *** ------------------------------------------------
  |											|
  |	Substitution edit:	Format: (reference/edit)	Example:(G/A)		|
  |	Insertion edit:		Format: (+insertion)		Example:(+ATCG)		|
  |	Deletion edit:		Format: (-deletion)		Example:(-ATCG)		|
  |											|
  ---------------------------------------------------------------------------------------

  *** Combination edit example *** ------------------------------------------------------
  |											|
  |	Reference:			ATGCTGTGAT G TCGTGATG    A			|
  |	Edit:				A--CTGTGAT C TCGTGATGatcgA			|
  |	Sequence format:	A(-TG)CTGTGAT(G/C)TCGTGATG(+atcg)A			|
  |											|
  ---------------------------------------------------------------------------------------[/blue]""")
        sys.exit()

    # Pilot_Screen():
    parser_pe_pilot_screen.add_argument("--pegRNAs", type=str, dest='pegRNAs_dir',help="Directory with pegRNAs from prime_designer() output", required=True)
    parser_pe_pilot_screen.add_argument("--mutations", type=str, dest='mutations_pt', help="Path to mutations file (COSMIC or ClinVar)", required=True)
    
    parser_pe_pilot_screen.add_argument("--database", type=str, choices=['COSMIC', 'ClinVar'], default='COSMIC', help="Database to use for priority mutations (Default: 'COSMIC')")
    parser_pe_pilot_screen.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    # epegRNA_linkers():
    parser_pe_epegRNA_linkers.add_argument('--pegRNAs', help='Path to pegRNAs file',required=True)

    parser_pe_epegRNA_linkers.add_argument('--epegRNA_motif_sequence', default='CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA', help='epegRNA motif sequence (Default: tevopreQ1)')
    parser_pe_epegRNA_linkers.add_argument('--linker_pattern', type=str, default=argparse.SUPPRESS, help='epegRNA linker pattern (Default: NNNNNNNN)')
    parser_pe_epegRNA_linkers.add_argument('--excluded_motifs', type=str, nargs="+", default=['Esp3I'], help="list of motifs or type IIS RE enzymes (i.e., Esp3I, BsaI, BspMI) to exclude from linker generation (Default: ['Esp3I'])")
    parser_pe_epegRNA_linkers.add_argument('--ckpt_dir', type=str, help='Checkpoint directory path (Default: ../epegRNAs/ckpt)', default='../epegRNAs/ckpt')
    parser_pe_epegRNA_linkers.add_argument('--ckpt_file', help='Checkpoint file name (Default: YYMMDD_HHMMSS_epegRNA_linkers.csv)', default=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_epegRNA_linkers.csv')
    parser_pe_epegRNA_linkers.add_argument('--ckpt_pt', type=str, default='', help='Previous checkpoint full path (Example: ../epegRNAs/ckpt/YYMMDD_HHMMSS_epegRNA_linkers.csv)')
    parser_pe_epegRNA_linkers.add_argument("--out_dir", type=str, help="Output directory (Default: ../epegRNAs)", default='../epegRNAs')
    parser_pe_epegRNA_linkers.add_argument("--out_file", type=str, help="Name of the output file (Default: epegRNAs.csv)", default='epegRNAs.csv')
    parser_pe_epegRNA_linkers.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    # merge():
    parser_pe_merge.add_argument("--epegRNAs", type=str, help="Directory or file with epegRNAs", required=True)
    parser_pe_merge.add_argument("--ngRNAs", type=str, help="Directory or file with ngRNAs", required=True)
    
    parser_pe_merge.add_argument("--ngRNAs_groups", type=str, dest='ngRNAs_groups_max', help="Maximum # of ngRNAs per epegRNA (Default: 3)", default=3)
    parser_pe_merge.add_argument("--epegRNA_suffix", type=str, help="Suffix for epegRNAs columns (Default: _epegRNA)", default='_epegRNA')
    parser_pe_merge.add_argument("--ngRNA_suffix", type=str, help="Suffix for ngRNAs columns (Default: _ngRNA)", default='_ngRNA')
    parser_pe_merge.add_argument("--out_dir", type=str, dest='dir', help="Output directory (Default: ../epeg_ngRNAs)", default='../epeg_ngRNAs')
    parser_pe_merge.add_argument("--out_file", type=str, dest='file', help="Name of the output file (Default: epeg_ngRNAs.csv)", default='epeg_ngRNAs.csv')
    parser_pe_merge.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    # sensor_designer():
    parser_pe_sensor_designer.add_argument("--pegRNAs", type=str, help="Path to pegRNAs file", required=True)
    
    parser_pe_sensor_designer.add_argument("--sensor_length", type=int, default=60, help="Total length of the sensor in bp (Default: 60)")
    parser_pe_sensor_designer.add_argument("--before_spacer", type=int, default=5, help="Amount of nucleotide context to put before the protospacer in the sensor (Default = 5)")
    parser_pe_sensor_designer.add_argument("--sensor_orientation", type=str, default='revcom', help="Orientation of the sensor relative to the protospacer (Options: 'revcom' [Default b/c minimize recombination] or ’forward’")
    parser_pe_sensor_designer.add_argument("--out_dir", type=str, help="Output directory (Default: ../pegRNAs_tester)", default='../sensor_designer')
    parser_pe_sensor_designer.add_argument("--out_file", type=str, help="Name of the output file (Default: pegRNAs.csv)", default='pegRNAs.csv')
    parser_pe_sensor_designer.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    # pegRNA_outcome():
    parser_pe_pegRNA_outcome.add_argument("--pegRNAs", type=str, help="Path to pegRNAs file", required=True)

    parser_pe_pegRNA_outcome.add_argument("--in_file", type=str, help="Path to PrimeDesign input file (required columns: target_name, target_sequence, index). Verify all expected edits are present in pegRNA library (for saturation mutagenesis only)")
    parser_pe_pegRNA_outcome.add_argument("--out_dir", type=str, help="Output directory (Default: ../pegRNA_outcome)", default='../pegRNA_outcome')
    parser_pe_pegRNA_outcome.add_argument("--out_file", type=str, help="Name of the output file (Default: pegRNAs.csv)", default='pegRNAs.csv')
    parser_pe_pegRNA_outcome.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    parser_pe_pegRNA_outcome.add_argument("--match_score", type=float, help="Match score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_outcome.add_argument("--mismatch_score", type=float, help="Mismatch score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_outcome.add_argument("--open_gap_score", type=float, help="Open gap score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_outcome.add_argument("--extend_gap_score", type=float, help="Extend gap score for pairwise alignment", default=argparse.SUPPRESS)
    
    # pegRNA_signature():
    parser_pe_pegRNA_signature.add_argument("--pegRNAs", type=str, help="Path to pegRNAs file", required=True)

    parser_pe_pegRNA_signature.add_argument("--config_key", type=str, help="Config file key (FWD primer_REV primer) with 'motif5' (flank5) & 'motif3' (flank3)", default=argparse.SUPPRESS)
    parser_pe_pegRNA_signature.add_argument("--flank5", type=str, dest='flank5_sequence', help=f"Flank 5' sequence", default=argparse.SUPPRESS)
    parser_pe_pegRNA_signature.add_argument("--flank3", type=str, dest='flank3_sequence', help=f"Flank 3' sequence", default=argparse.SUPPRESS)
    parser_pe_pegRNA_signature.add_argument("--flank5_length", type=int, help="Length of flank5 sequence to include in alignment if provided (Default: 0)", default=0)
    parser_pe_pegRNA_signature.add_argument("--flank3_length", type=int, help="Length of flank3 sequence to include in alignment if provided (Default: 0)", default=0)
    parser_pe_pegRNA_signature.add_argument("--reference_sequence", type=str, help="Column name for reference sequences (Default: 'Reference_sequence')", default='Reference_sequence')
    parser_pe_pegRNA_signature.add_argument("--edit_sequence", type=str, help="Column name for edit sequences (Default: 'Edit_sequence')", default='Edit_sequence')
    parser_pe_pegRNA_signature.add_argument("--out_dir", type=str, help="Output directory (Default: ../pegRNA_signature)", default='../pegRNA_signature')
    parser_pe_pegRNA_signature.add_argument("--out_file", type=str, help="Name of the output file (Default: pegRNAs.csv)", default='pegRNAs.csv')
    parser_pe_pegRNA_signature.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    parser_pe_pegRNA_signature.add_argument("--match_score", type=float, help="Match score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_signature.add_argument("--mismatch_score", type=float, help="Mismatch score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_signature.add_argument("--open_gap_score", type=float, help="Open gap score for pairwise alignment", default=argparse.SUPPRESS)
    parser_pe_pegRNA_signature.add_argument("--extend_gap_score", type=float, help="Extend gap score for pairwise alignment", default=argparse.SUPPRESS)

    # epegRNA_fastas():
    parser_pe_epegRNA_fasta.add_argument("--df", type=str, help="Path to epegRNAs file", required=True)
    parser_pe_epegRNA_fasta.add_argument("--linearized_vector", type=str, help="Linearized vector sequence such that final plasmid is linearized vector + insert (string or .fasta)", required=True)

    parser_pe_epegRNA_fasta.add_argument("--out_dir", type=str, help="Output directory for FASTA files (Default: ./fasta)", default='./fasta')
    parser_pe_epegRNA_fasta.add_argument("--id", type=str, help="epegRNA ID column name (Default: ID)", default='ID')
    parser_pe_epegRNA_fasta.add_argument("--dont_tG", dest="tG", default=True, action="store_false", help="Don't add 5' G to spacer if needed")
    parser_pe_epegRNA_fasta.add_argument("--dont_make_extension", dest="make_extension", default=True, action="store_false", help="Don't build extension from RTT, PBS, and linker")
    parser_pe_epegRNA_fasta.add_argument("--epegRNA_spacer", type=str, help="epegRNA spacer column name (Default: Spacer_sequence)", default='Spacer_sequence')
    parser_pe_epegRNA_fasta.add_argument("--epegRNA_scaffold", type=str, default="Scaffold_sequence", help="epegRNA scaffold column")
    parser_pe_epegRNA_fasta.add_argument("--epegRNA_extension", type=str, default="Extension_sequence", help="epegRNA extension column")
    parser_pe_epegRNA_fasta.add_argument("--epegRNA_RTT", type=str, default="RTT_sequence", help="epegRNA RTT column name")
    parser_pe_epegRNA_fasta.add_argument("--epegRNA_PBS", type=str, default="PBS_sequence", help="epegRNA PBS column name")
    parser_pe_epegRNA_fasta.add_argument("--epegRNA_linker", type=str, default="Linker_sequence", help="epegRNA Linker column name")
    parser_pe_epegRNA_fasta.add_argument("--no_literals", action='store_false', dest='literal_eval', help="Do not convert string representations", default=True)

    # Set defaults
    parser_pe_prime_designer.set_defaults(func=pe.prime_designer)
    parser_pe_pilot_screen.set_defaults(func=pe.pilot_screen)
    parser_pe_epegRNA_linkers.set_defaults(func=pe.epegRNA_linkers)
    parser_pe_merge.set_defaults(func=pe.merge)
    parser_pe_sensor_designer.set_defaults(func=pe.sensor_designer)
    parser_pe_pegRNA_outcome.set_defaults(func=pe.pegRNA_outcome)
    parser_pe_pegRNA_signature.set_defaults(func=pe.pegRNA_signature)
    parser_pe_epegRNA_fasta.set_defaults(func=pe.epegRNA_fasta)