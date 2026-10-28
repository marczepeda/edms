# Design pegRNAs and ngRNAs for prime editing
##### Import libraries
import os
import sys
import re
import difflib
import functools
import time
import argparse
import logging
from argparse import RawTextHelpFormatter


logger = logging.getLogger(__name__)

# Module-level placeholders populated by _configure_cli(); only meaningful after main() runs.
file_in = None
pe_format = None
pbs_length_list = 0
rtt_length_list = 0
nicking_distance_minimum = 0
nicking_distance_maximum = 120
filter_c1_extension = False
filter_homopolymer_ts = False
silent_mutation = False
genome_wide_design = False
saturation_mutagenesis = False
number_of_pegrnas = 0
number_of_silent_mutations = 0
number_of_ngrnas = 1
nicking_distance_pooled = 75
homology_downstream = 15
pbs_length_pooled = 14
rtt_max_length_pooled = 50
out_dir = None


def _configure_cli():
	"""Parse CLI args and configure logging/output dir. Only called from main(), never on import."""
	global args, file_in, pe_format, pbs_length_list, rtt_length_list, nicking_distance_minimum, nicking_distance_maximum
	global filter_c1_extension, filter_homopolymer_ts, silent_mutation, genome_wide_design, saturation_mutagenesis
	global number_of_pegrnas, number_of_silent_mutations, number_of_ngrnas, nicking_distance_pooled, homology_downstream, pbs_length_pooled, rtt_max_length_pooled
	global out_dir, logger

	parser = argparse.ArgumentParser(description = '''----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	Software for the design of pegRNAs for flexible prime editing! Please visit ----- https://github.com/jyhsu15/PrimeDesign ----- for more documentation on how to use the software.\nModified by Marc Zepeda for integration into the EDMS pipeline. Please visit ----- https://github.com/marczepeda/edms ----- for more documentation on how to use the software.
	----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------''', formatter_class=RawTextHelpFormatter)

	# Inputs for de-novo design of pegRNAs and nicking gRNAs
	parser.add_argument('-f', '--file', required = True, type = str, help = '''Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (Required)

	*** Example .TXT file *** --------------------------------------------------------------
	|											|
	|	target_01_substitution	ATGTGCTGTGATGGTAT(G/A)CCGGCGTAGTAATCGTAGC		|
	|	target_01_insertion	ATGTGCTGTGATGGTATG(+ATCTCGATGA)CCGGCGTAGTAATCGTAGC	|
	|	target_01_deletion	ATGTGCTGTGATGG(-TATGCCG)GCGTAGTAATCGTAGC		|
	|											|
	 ---------------------------------------------------------------------------------------

	*** Example .CSV file *** --------------------------------------------------------------
	|											|
	|	target_01_substitution,ATGTGCTGTGATGGTAT(G/A)CCGGCGTAGTAATCGTAGC		|
	|	target_01_insertion,ATGTGCTGTGATGGTATG(+ATCTCGATGA)CCGGCGTAGTAATCGTAGC		|
	|	target_01_deletion,ATGTGCTGTGATGG(-TATGCCG)GCGTAGTAATCGTAGC			|
	|											|
	 ---------------------------------------------------------------------------------------

	*** Formatting different DNA edits *** -------------------------------------------------
	|											|
	|	Substitution edit:	Format: (reference/edit)	Example:(G/A)		|
	|	Insertion edit:		Format: (+insertion)		Example:(+ATCG)		|
	|	Deletion edit:		Format: (-deletion)		Example:(-ATCG)		|
	|											|
	 ---------------------------------------------------------------------------------------

	*** Combination edit example *** -------------------------------------------------------
	|											|
	|	Reference:			ATGCTGTGAT G TCGTGATG    A			|
	|	Edit:				A--CTGTGAT C TCGTGATGatcgA			|
	|	Sequence format:	A(-TG)CTGTGAT(G/C)TCGTGATG(+atcg)A			|
	|											|
	 ---------------------------------------------------------------------------------------

	''')

	# Inputs for the design parameters of pegRNAs and nicking gRNAs
	parser.add_argument('-pe_format', '--pe_format', type = str, default = 'NNNNNNNNNNNNNNNNN/NNN[NGG]', help = "***** Prime editing formatting including the spacer, cut index -> /, and protospacer adjacent motif (PAM) -> [PAM] (Default: NNNNNNNNNNNNNNNNN/NNN[NGG]). Examples: NNNNNNNNNNNNNNNNN/NNN[NGG], NNNNNNNNNNNNNNNNN/NNN[NG] *****\n\n")
	parser.add_argument('-pbs', '--pbs_length_list', type = int, default = 0, nargs = '+', help = '***** List of primer binding site (PBS) lengths for the pegRNA extension (Default: 10 to 16 nt). Example: 12 13 14 15 *****\n\n')
	parser.add_argument('-rtt', '--rtt_length_list', type = int, default = 0, nargs = '+', help = '***** List of reverse transcription (RT) template lengths for the pegRNA extension (Default: 10 to 50 nt). Example: 10 15 20 *****\n')
	parser.add_argument('-nick_dist_min', '--nicking_distance_minimum', type = int, default = 0, nargs = '+', help = '***** Minimum nicking distance for designing ngRNAs upstream and downstream of a pegRNA (Default: 0). *****\n\n')
	parser.add_argument('-nick_dist_max', '--nicking_distance_maximum', type = int, default = 120, nargs = '+', help = '***** Maximum nicking distance for designing ngRNAs upstream and downstream of a pegRNA (Default: 100). *****\n\n')
	parser.add_argument('-filter_c1', '--filter_c1_extension', action='store_true', help = '***** Option to filter against pegRNA extensions that start with a C base. *****\n\n')
	parser.add_argument('-filter_homopolymer_ts', '--filter_homopolymer_ts', action='store_true', help = '***** Option to filter out spacer sequences with homopolymer Ts (>3). *****\n\n')
	parser.add_argument('-silent_mut', '--silent_mutation', action='store_true', help = '***** Introduce silent mutation into or around the PAM assuming the sequence is in-frame. Currently only available with SpCas9 PE (i.e., pe_format = NNNNNNNNNNNNNNNNN/NNN[NGG]). *****\n\n')
	parser.add_argument('-genome_wide', '--genome_wide_design', action='store_true', help = '***** Whether or not this is a genome-wide pooled design. This option designs a set of pegRNAs per input without ranging PBS and RTT parameters. (Default: False) *****\n\n')
	parser.add_argument('-sat_mut', '--saturation_mutagenesis', default = False, choices = ['aa', 'aa_subs', 'aa_ins', 'aa_dels', 'aa_silent', 'base'], type = str, help = '***** Saturation mutagenesis design with prime editing. The \'aa\' option makes all amino acid substitutions (\'aa_subs\'),  +1 amino acid insertions (\'aa_ins\'), and -1 amino acid deletions (\'aa_dels\'). The \'aa_silent\' option makes all possible silent codon changes. The \'base\' option makes DNA base changes. (Default: False) *****\n\n')
	parser.add_argument('-n_pegrnas', '--number_of_pegrnas', default = 0, type = int, help = '***** The maximum number of pegRNAs to design for each input sequence. 0 keeps every pegRNA meeting the criteria. The pegRNAs are ranked by 1) PAM disrupted > Spacer disrupted > PAM & Spacer intact then 2) distance to edit. (Default: 0) *****\n\n')
	parser.add_argument('-n_silent_mut', '--number_of_silent_mutations', default = 0, type = int, help = '***** The maximum number of silent mutation options to design for each pegRNA when the input file supplies the edits (not saturation mutagenesis). 0, the default, keeps a curated set: every option disrupting the PAM or the spacer, plus the non-disrupting option closest to the PAM. A non-zero value instead takes that many from the full ranked list. Every pegRNA additionally gets a design with no silent mutation. (Default: 0) *****\n\n')
	parser.add_argument('-n_ngrnas', '--number_of_ngrnas', default = 1, type = int, help = '***** The maximum number of ngRNAs to design for each input sequence. The ngRNAs are ranked by 1) PE3b-seed > PE3b-nonseed > PE3 then 2) deviation from nicking_distance_pooled. (Default: 1) *****\n\n')
	parser.add_argument('-nick_dist_pooled', '--nicking_distance_pooled', default = 75, type = int, help = '***** The nicking distance between pegRNAs and ngRNAs for pooled designs. PE3b annotation is priority, followed by nicking distance closest to this parameter. (Default: 75 bp) *****\n\n')
	parser.add_argument('-homology_downstream', '--homology_downstream', default = 15, type = int, help = '***** This parameter determines the minimum RT extension length downstream of an edit for pegRNA designs. (Default: 15) *****\n\n')
	parser.add_argument('-pbs_pooled', '--pbs_length_pooled', type = int, default = 14, help = '***** The PBS length to design pegRNAs for pooled design applications. (Default: 14 nt) *****\n\n')
	parser.add_argument('-rtt_pooled', '--rtt_max_length_pooled', type = int, default = 50, help = '***** The maximum RTT length to design pegRNAs for pooled design applications. (Default: 50 nt) *****\n\n')


	# Output directory
	parser.add_argument('-out', '--out_dir', default = './PrimeDesign/DATETIMESTAMP_PrimeDesign', type = str, help = '***** Name of output directory (Default: ./DATETIMESTAMP_PrimeDesign). *****\n\n')

	args = parser.parse_args()

	##### Initialize arguments
	file_in = args.file

	pe_format = args.pe_format
	pbs_length_list = args.pbs_length_list
	rtt_length_list = args.rtt_length_list
	nicking_distance_minimum = args.nicking_distance_minimum
	nicking_distance_maximum = args.nicking_distance_maximum
	filter_c1_extension = args.filter_c1_extension
	filter_homopolymer_ts = args.filter_homopolymer_ts
	silent_mutation = args.silent_mutation

	genome_wide_design = args.genome_wide_design
	saturation_mutagenesis = args.saturation_mutagenesis
	number_of_pegrnas = args.number_of_pegrnas
	number_of_silent_mutations = args.number_of_silent_mutations
	number_of_ngrnas = args.number_of_ngrnas
	nicking_distance_pooled = args.nicking_distance_pooled

	homology_downstream = args.homology_downstream
	pbs_length_pooled = args.pbs_length_pooled
	rtt_max_length_pooled = args.rtt_max_length_pooled

	# Default PBS and RTT lengths to design
	if pbs_length_list == 0:
		pbs_length_list = list(range(10, 16))
	if rtt_length_list == 0:
		rtt_length_list = list(range(10, 50))

	# Output directory date and time stamped
	out_dir = args.out_dir
	if out_dir == './PrimeDesign/DATETIMESTAMP_PrimeDesign':
		out_dir = './PrimeDesign/%s_PrimeDesign' % str(time.strftime("%y%m%d_%H.%M.%S", time.localtime()))

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	# Initialize logger
	logger = logging.getLogger()
	logger.setLevel(logging.DEBUG)

	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

	fh = logging.FileHandler(out_dir + '/PrimeDesign.log')
	fh.setLevel(logging.DEBUG)
	fh.setFormatter(formatter)
	logger.addHandler(fh)

	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	ch.setFormatter(formatter)
	logger.addHandler(ch)
##### IUPAC code map
iupac2bases_dict = {'A':'A','T':'T','C':'C','G':'G','a':'a','t':'t','c':'c','g':'g',
'R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]','K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACTG]',
'r':'[ag]','y':'[ct]','s':'[gc]','w':'[at]','k':'[gt]','m':'[ac]','b':'[cgt]','d':'[agt]','h':'[act]','v':'[acg]','n':'[actg]',
'(':'(',')':')','+':'+','-':'-','/':'/'}

def iupac2bases(iupac):

	try:
		bases = iupac2bases_dict[iupac]
	except:
		logger.error('Symbol %s is not within the IUPAC nucleotide code ...' % str(iupac))
		sys.exit(1)

	return(bases)

# GC content
def gc_content(sequence):
	sequence = sequence.upper()
	GC_count = sequence.count('G') + sequence.count('C')
	GC_content = float(GC_count)/float(len(sequence))

	return("%.2f" % GC_content)

# Reverse complement function
def reverse_complement(sequence):
	sequence = sequence
	new_sequence = ''
	for base in sequence:
		if base == 'A':
			new_sequence += 'T'
		elif base == 'T':
			new_sequence += 'A'
		elif base == 'C':
			new_sequence += 'G'
		elif base == 'G':
			new_sequence += 'C'
		elif base == 'a':
			new_sequence += 't'
		elif base == 't':
			new_sequence += 'a'
		elif base == 'c':
			new_sequence += 'g'
		elif base == 'g':
			new_sequence += 'c'
		elif base == '[':
			new_sequence += ']'
		elif base == ']':
			new_sequence += '['
		elif base == '+':
			new_sequence += '+'
		elif base == '-':
			new_sequence += '-'
		elif base == '/':
			new_sequence += '/'
		elif base == '(':
			new_sequence += ')'
		elif base == ')':
			new_sequence += '('
	return(new_sequence[::-1])

# Amino acid code
codon_dict = {
	'GGG':['Gly','G', 0.25],'GGA':['Gly','G', 0.25],'GGT':['Gly','G', 0.16],'GGC':['Gly','G', 0.34],
	'GAG':['Glu','E', 0.58],'GAA':['Glu','E', 0.42],'GAT':['Asp','D', 0.46],'GAC':['Asp','D', 0.54],
	'GTG':['Val','V', 0.47],'GTA':['Val','V', 0.11],'GTT':['Val','V', 0.18],'GTC':['Val','V', 0.24],
	'GCG':['Ala','A', 0.11],'GCA':['Ala','A', 0.23],'GCT':['Ala','A', 0.26],'GCC':['Ala','A', 0.4],
	'AGG':['Arg','R', 0.2],'AGA':['Arg','R', 0.2],'AGT':['Ser','S', 0.15],'AGC':['Ser','S', 0.24],
	'AAG':['Lys','K', 0.58],'AAA':['Lys','K', 0.42],'AAT':['Asn','N', 0.46],'AAC':['Asn','N', 0.54],
	'ATG':['Met','M', 1],'ATA':['Ile','I', 0.16],'ATT':['Ile','I', 0.36],'ATC':['Ile','I', 0.48],
	'ACG':['Thr','T', 0.12],'ACA':['Thr','T', 0.28],'ACT':['Thr','T', 0.24],'ACC':['Thr','T', 0.36],
	'TGG':['Trp','W', 1],'TGA':['End','X', 0.52],'TGT':['Cys','C', 0.45],'TGC':['Cys','C', 0.55],
	'TAG':['End','X', 0.2],'TAA':['End','X', 0.28],'TAT':['Tyr','Y', 0.43],'TAC':['Tyr','Y', 0.57],
	'TTG':['Leu','L', 0.13],'TTA':['Leu','L', 0.07],'TTT':['Phe','F', 0.45],'TTC':['Phe','F', 0.55],
	'TCG':['Ser','S', 0.06],'TCA':['Ser','S', 0.15],'TCT':['Ser','S', 0.18],'TCC':['Ser','S', 0.22],
	'CGG':['Arg','R', 0.21],'CGA':['Arg','R', 0.11],'CGT':['Arg','R', 0.08],'CGC':['Arg','R', 0.19],
	'CAG':['Gln','Q', 0.75],'CAA':['Gln','Q', 0.25],'CAT':['His','H', 0.41],'CAC':['His','H', 0.59],
	'CTG':['Leu','L', 0.41],'CTA':['Leu','L', 0.07],'CTT':['Leu','L', 0.13],'CTC':['Leu','L', 0.2],
	'CCG':['Pro','P', 0.11],'CCA':['Pro','P', 0.27],'CCT':['Pro','P', 0.28],'CCC':['Pro','P', 0.33],
}

# Create codon swap dictionaries
aa2codon = {}
for codon in codon_dict:
	if codon_dict[codon][1] not in aa2codon:
		aa2codon[codon_dict[codon][1]] = []

	aa2codon[codon_dict[codon][1]].append([codon, codon_dict[codon][2]])

for codon in aa2codon:
	aa2codon[codon] = sorted(aa2codon[codon], key = lambda x: x[1], reverse = True)

def codon_nt_difference(codon_ref, codon_edit):
	return(sum([1 if ref_base != edit_base else 0 for ref_base, edit_base in zip(codon_ref.upper(), codon_edit.upper())]))

def silent_codon_in_allowed_range(target_design_entry, codon_start_abs, codon_end_abs):
	"""Return True when a candidate silent codon is inside the saturation-mutagenesis region.

	The range is only defined for saturation mutagenesis. For normal PrimeDesign
	inputs, this helper returns True and does not change behavior.
	Coordinates are 0-indexed and half-open: [start, end).
	"""
	allowed_range = target_design_entry.get('silent_mutation_allowed_range', None)
	if allowed_range is None:
		return True
	return codon_start_abs >= allowed_range[0] and codon_end_abs <= allowed_range[1]

def silent_mutation_edit_sequence_in_allowed_range(target_design_entry, edit_sequence, silent_mutation_edit_sequence):
	"""Validate that any additional silent codon changes remain inside the allowed region.

	This is a final safety gate for saturation mutagenesis. It compares the edit
	sequence before and after adding the silent PAM/synonymous mutation. Any
	new base differences are assigned to their codon and the whole codon must be
	inside silent_mutation_allowed_range.
	"""
	allowed_range = target_design_entry.get('silent_mutation_allowed_range', None)
	if allowed_range is None or silent_mutation_edit_sequence in ['', None]:
		return True
	if len(edit_sequence) != len(silent_mutation_edit_sequence):
		return False
	changed_codons = set()
	for i, (base_before, base_after) in enumerate(zip(edit_sequence.upper(), silent_mutation_edit_sequence.upper())):
		if base_before != base_after:
			codon_start = int(i / 3) * 3
			changed_codons.add((codon_start, codon_start + 3))
	return all(silent_codon_in_allowed_range(target_design_entry, start, end) for start, end in changed_codons)

# Written to Silent_mutation_relative_to_edit when a pegRNA carries no silent mutation.
# Deliberately not the string 'None': that is in pandas' default NA list, so read_csv turns
# it into NaN and the row looks as though the value were missing rather than meaningful.
NO_SILENT_MUTATION = 'no_silent_mutation'

def ranked_silent_mutations(target_design_entry, edit_sequence, pe_nick_edit_idx, rtt_length,
	strand, pam_search, edit_start_abs, edit_stop_abs, homology_downstream = 0):
	"""Enumerate every silent codon substitution the RTT can encode, best candidate first.

	Replaces the previous approach of inspecting only the one or two codons adjacent to
	the PAM: every in-frame codon lying wholly inside the RTT (and, under saturation
	mutagenesis, inside the allowed region) is considered.

	Candidates are ranked by, in order:
		1. disrupts the PAM
		2. disrupts the portion of the spacer the RTT can reach
		3. number of nucleotide differences from the reference codon, most first
		4. human codon usage of the replacement codon, most frequent first

	The last homology_downstream nucleotides of the RTT are excluded. That stretch is the 3'
	homology arm the reverse-transcribed flap anneals through, so a mismatch there costs
	incorporation efficiency; it is the one part of the RTT that must match the genome.

	Coordinates are 0-indexed, half-open, and in edit_sequence space. The RTT runs 3' of
	the nick on the (+) strand and 5' of the nick on the (-) strand, so the PAM and the
	reachable spacer bases sit at mirrored offsets.
	"""
	if strand == '+':
		rtt_lo, rtt_hi = pe_nick_edit_idx, pe_nick_edit_idx + rtt_length - homology_downstream
		pam_lo, pam_hi = pe_nick_edit_idx + 3, pe_nick_edit_idx + 6
		spacer_lo, spacer_hi = pe_nick_edit_idx, pe_nick_edit_idx + 3
		pam_pattern = pam_search
	else:
		rtt_lo, rtt_hi = pe_nick_edit_idx - rtt_length + homology_downstream, pe_nick_edit_idx
		pam_lo, pam_hi = pe_nick_edit_idx - 6, pe_nick_edit_idx - 3
		spacer_lo, spacer_hi = pe_nick_edit_idx - 3, pe_nick_edit_idx
		pam_pattern = reverse_complement(pam_search)

	rtt_lo = max(rtt_lo, 0)
	rtt_hi = min(rtt_hi, len(edit_sequence))
	pam_intact_before = bool(re.search(pam_pattern, edit_sequence[pam_lo:pam_hi].upper())) if pam_lo >= 0 else False

	candidates = []
	first_codon = rtt_lo + ((3 - rtt_lo % 3) % 3)  # First in-frame codon start at or after rtt_lo
	for codon_start in range(first_codon, rtt_hi - 2, 3):
		codon_end = codon_start + 3

		if codon_start < edit_stop_abs and codon_end > edit_start_abs:
			continue  # Codon carries part of the intended edit
		if not silent_codon_in_allowed_range(target_design_entry, codon_start, codon_end):
			continue  # Outside the saturation mutagenesis region

		codon_ref = edit_sequence[codon_start:codon_end].upper()
		if codon_ref not in codon_dict:
			continue  # Ambiguous or incomplete codon
		aa_identity = codon_dict[codon_ref][1]

		for codon_substitute in aa2codon[aa_identity]:
			codon_edit = codon_substitute[0].upper()
			if codon_edit == codon_ref:
				continue

			new_edit_sequence = edit_sequence[:codon_start] + codon_edit.lower() + edit_sequence[codon_end:]
			disrupts_pam = pam_intact_before and not re.search(pam_pattern, new_edit_sequence[pam_lo:pam_hi].upper())
			disrupts_spacer = any(edit_sequence[i].upper() != new_edit_sequence[i].upper()
				for i in range(max(spacer_lo, codon_start), min(spacer_hi, codon_end)))

			# 0 when the codon overlaps the PAM, otherwise the gap between them
			distance_to_pam = max(0, pam_lo - codon_end, codon_start - pam_hi)

			candidates.append({
				'distance_to_pam': distance_to_pam,
				'codon_start_abs': codon_start,
				'codon_end_abs': codon_end,
				'aa_number': int(codon_start / 3) + 1,
				'codon_ref': codon_ref,
				'codon_edit': codon_edit,
				'aa_identity': aa_identity,
				'usage_ref': codon_dict[codon_ref][2],
				'usage_edit': codon_dict[codon_edit][2],
				'nt_difference': codon_nt_difference(codon_ref, codon_edit),
				'disrupts_pam': bool(disrupts_pam),
				'disrupts_spacer': bool(disrupts_spacer),
				'silent_mutation_edit_sequence': new_edit_sequence,
			})

	return sorted(candidates, key = lambda c: (-c['disrupts_pam'], -c['disrupts_spacer'],
		-c['nt_difference'], -c['usage_edit'], c['codon_edit']))

SILENT_CODON_COLUMNS = [
	'Target_reference_codon', 'Target_reference_AA', 'Target_reference_codon_usage',
	'Target_reference_AA_number', 'Target_reference_codon_index',
	'Target_edit_codon', 'Target_edit_AA', 'Target_edit_codon_usage',
	'Target_edit_AA_number', 'Target_edit_codon_index',
	'Silent_reference_codon', 'Silent_reference_AA', 'Silent_reference_codon_usage',
	'Silent_reference_AA_number', 'Silent_reference_codon_index',
	'Silent_edit_codon', 'Silent_edit_AA', 'Silent_edit_codon_usage',
	'Silent_edit_AA_number', 'Silent_edit_codon_index',
]

def target_aa_number(target_design_entry, codon_start):
	"""Amino acid number for an in-frame codon, matching the Edit column convention.

	Numbering runs from the user-supplied index in the third input column. Under saturation
	mutagenesis that index names the first residue of the parenthesized region; otherwise it
	names the first residue of the target sequence. With no index the numbering starts at 1.
	"""
	allowed_range = target_design_entry.get('silent_mutation_allowed_range', None)
	origin = allowed_range[0] if allowed_range is not None else 0
	return (codon_start - origin) // 3 + int(target_design_entry.get('aa_index', 1))

def codon_annotation(sequence, codon_start, aa_number):
	"""[codon, amino acid, human codon usage, AA number, nt index] for one in-frame codon.

	nt index is the 0-indexed start of the codon in the target sequence, so it indexes
	Reference_sequence and Edit_sequence directly.
	"""
	codon = str(sequence)[codon_start:codon_start + 3].upper()
	if len(codon) != 3 or codon not in codon_dict:
		return ['', '', '', '', '']
	return [codon, codon_dict[codon][1], codon_dict[codon][2], aa_number, codon_start]

def codon_columns(target_design_entry, silent_selected):
	"""The 20 codon columns: the intended edit's codon, then the silent mutation's codon.

	Each is reported before and after the change as
	[nucleotides, amino acid, codon usage, AA number, nt index]. Blank when the codon is
	not resolvable (indel edits shift the frame) or when no silent mutation was made.
	"""
	edit_start = target_design_entry.get('edit_start_in_ref', None)
	if edit_start is None:
		target_cols = [''] * 10
	else:
		codon_start = int(int(edit_start) / 3) * 3
		number = target_aa_number(target_design_entry, codon_start)
		target_cols = (codon_annotation(target_design_entry['reference_sequence'], codon_start, number) +
			codon_annotation(target_design_entry['edit_sequence'], codon_start, number))

	if not silent_selected:
		return target_cols + [''] * 10

	start = silent_selected['codon_start_abs']
	number = target_aa_number(target_design_entry, start)
	ref, edit = silent_selected['codon_ref'], silent_selected['codon_edit']
	return target_cols + [
		ref, codon_dict[ref][1], codon_dict[ref][2], number, start,
		edit, codon_dict[edit][1], codon_dict[edit][2], number, start]

def extend_rtt_for_silent_mutations(target_design_entry, edit_sequence, pe_nick_edit_idx, rtt_length,
	strand, pam_search, edit_start_abs, edit_stop_abs, homology_downstream, rtt_max_length):
	"""Grow the RTT until it can carry a silent mutation outside the 3' homology arm.

	Saturation mutagenesis sizes the RTT as the distance to the end of the edit plus the
	homology arm, so once the arm is excluded there is nothing left between the two whenever
	the nick sits close to the edit. Lengthening the RTT pushes the arm 3' and opens more
	positions, so it is grown one nucleotide at a time until a candidate appears or the length
	would reach rtt_max_length. A single nucleotide is enough to matter when the window
	boundary falls mid-codon, and stepping one at a time keeps the extension as short as
	possible, which a codon-sized step would not.

	Returns (rtt_length, candidates). The original length comes back unchanged when no
	extension is needed, and also when no extension produces a candidate, so a pegRNA is
	never lengthened for nothing.
	"""
	def candidates_for(length):
		return ranked_silent_mutations(target_design_entry, edit_sequence, pe_nick_edit_idx,
			length, strand, pam_search, edit_start_abs, edit_stop_abs, homology_downstream)

	candidates = candidates_for(rtt_length)
	if candidates:
		return rtt_length, candidates

	extended = rtt_length + 1
	while extended < rtt_max_length:
		candidates = candidates_for(extended)
		if candidates:
			return extended, candidates
		extended += 1
	return rtt_length, []

def build_ngRNA_entries(target_design_entry, ngRNA_strand, allele, pe_nick_ref_idx, silent_mutation,
	pe_format, nicking_distance_minimum, nicking_distance_maximum, nicking_distance_pooled,
	ranked = True):
	"""Ranked ngRNA designs for one pegRNA allele.

	PE3b ngRNAs bind the edited strand, so their spacer has to be cut from the allele the
	pegRNA actually installs. Building them per allele stops a pegRNA carrying no silent
	mutation from inheriting spacers computed for one that does, which would not bind.

	ngRNA_strand is the strand the ngRNAs target, opposite the pegRNA's. The seed is the
	PAM-proximal ten nucleotides, which sit at the start of a (-) strand spacer and the end
	of a (+) strand one, and (-) strand spacers are stored reverse complemented.

	ranked matches the two entry shapes already in use: pooled designs lead each entry with a
	sort key combining the PE3b annotation and the deviation from nicking_distance_pooled, and
	are returned in that order; regular designs carry no key and keep the search order.
	"""
	entries = []
	for ngRNA in target_design_entry['ngRNA'][ngRNA_strand]:
		(ng_nick_ref_idx, ng_edit_start_idx, ng_edit_end_idx, ng_full_search_edit,
			ng_spacer_sequence_edit, ng_pam_edit, ng_annotate) = ngRNA
		nick_distance = ng_nick_ref_idx - pe_nick_ref_idx

		if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]') and len(allele) > 0:
			ng_spacer_sequence_edit = allele[ng_edit_start_idx:ng_edit_end_idx]

			# lowercase marks the edited bases
			mutation_indices = [i for i, base in enumerate(ng_spacer_sequence_edit) if base.islower()]
			if mutation_indices:
				in_seed = [i for i in mutation_indices if (i < 10 if ngRNA_strand == '-' else i >= 10)]
				ng_annotate = 'PE3b-seed' if in_seed else 'PE3b-nonseed'
			else:
				ng_annotate = 'PE3'

		if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):
			if ngRNA_strand == '-':
				entry = [ng_nick_ref_idx, reverse_complement(ng_spacer_sequence_edit),
					reverse_complement(ng_pam_edit), ng_annotate, '-', nick_distance]
			else:
				entry = [ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit,
					ng_annotate, '+', nick_distance]

			if ranked:
				ng_code = {'PE3b-seed': 0, 'PE3b-nonseed': 1}.get(ng_annotate, 2)
				offset = abs(abs(nick_distance) - nicking_distance_pooled)
				entry = [str(ng_code) + '0' * (3 - len(str(offset))) + str(offset)] + entry
			entries.append(entry)

	return sorted(entries) if ranked else entries

def select_silent_mutations(silent_candidates, maximum = 0):
	"""Which silent mutations a pegRNA should encode.

	With maximum 0, the default, a curated set is returned: every candidate that disrupts the
	PAM or the spacer, plus the one non-disrupting candidate closest to the PAM. Disruptors
	block re-cutting, and the closest non-disrupting codon is the most useful of the rest,
	since a change near the PAM is the most likely to interfere with re-binding.

	A non-zero maximum instead takes the top N of the full ranked list, which is how to ask
	for more options than the curated set, non-disrupting ones included.

	[None] means emit one pegRNA with no silent mutation at all.
	"""
	if maximum:
		return list(silent_candidates[:maximum]) or [None]

	disrupting = [candidate for candidate in silent_candidates
		if candidate['disrupts_pam'] or candidate['disrupts_spacer']]
	other = [candidate for candidate in silent_candidates
		if not (candidate['disrupts_pam'] or candidate['disrupts_spacer'])]

	# min() keeps the first of any tie, which is the better one by the existing ranking
	closest = [min(other, key = lambda candidate: candidate['distance_to_pam'])] if other else []
	return disrupting + closest or [None]

def with_no_silent_mutation_option(silent_options):
	"""Append the no-silent-mutation design so it is offered alongside the silent ones.

	A pegRNA carrying no silent mutation is still worth having, so every design emits one in
	addition to its silent options. It is written with Silent_mutation_relative_to_edit as
	'None'. Returns the list unchanged when it already is the no-silent-mutation case.
	"""
	return list(silent_options) if silent_options == [None] else list(silent_options) + [None]

def pam_spacer_annotation(reference_sequence, allele, pe_nick_ref_idx, pe_nick_edit_idx, strand,
	pam_search, cut_idx, spacer_end_idx, pam_start_idx, pam_end_idx):
	"""Classify what the pegRNA can no longer re-cut in the final allele.

		'PAM_disrupted'        the PAM no longer matches the PAM pattern
		'Spacer_disrupted'     the PAM still matches, but spacer bases 3' of the cut have changed
		'PAM_&_Spacer_Intact'  neither

	The protospacer core 5' of the cut must match the edited allele for the pegRNA to exist
	at all, so the only spacer bases that can differ are those between the cut and the PAM.
	Reference and allele are each read around their own nick in (+) orientation, which keeps
	the comparison valid when an indel shifts downstream coordinates.
	"""
	pam_lo, pam_hi = pam_start_idx - cut_idx, pam_end_idx - cut_idx
	tail = spacer_end_idx - cut_idx
	if strand == '+':
		allele_pam = allele[pe_nick_edit_idx + pam_lo:pe_nick_edit_idx + pam_hi]
		allele_tail = allele[pe_nick_edit_idx:pe_nick_edit_idx + tail]
		reference_tail = reference_sequence[pe_nick_ref_idx:pe_nick_ref_idx + tail]
		pattern = pam_search
	else:
		allele_pam = allele[pe_nick_edit_idx - pam_hi:pe_nick_edit_idx - pam_lo]
		allele_tail = allele[pe_nick_edit_idx - tail:pe_nick_edit_idx]
		reference_tail = reference_sequence[pe_nick_ref_idx - tail:pe_nick_ref_idx]
		pattern = reverse_complement(pam_search)

	if not re.search(pattern, allele_pam.upper()):
		return 'PAM_disrupted'
	if allele_tail.upper() != reference_tail.upper():
		return 'Spacer_disrupted'
	return 'PAM_&_Spacer_Intact'

NUCLEOTIDE_CHANGE_COLUMNS = ['Edit_sequence_change_indexes', 'Reference_sequence_change_indexes']

@functools.lru_cache(maxsize=8192)
def _nucleotide_change_indexes(reference_sequence, edit_sequence):
	"""Indexes of every nucleotide difference, in edit and in reference coordinates."""
	reference_sequence = reference_sequence.upper()
	edit_sequence = edit_sequence.upper()

	# Equal lengths means substitutions only, so the two coordinate systems coincide
	if len(reference_sequence) == len(edit_sequence):
		changed = tuple(i for i, (a, b) in enumerate(zip(reference_sequence, edit_sequence)) if a != b)
		return changed, changed

	edit_positions = []
	reference_positions = []
	for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(a=reference_sequence, b=edit_sequence, autojunk=False).get_opcodes():
		if tag == 'equal':
			continue
		reference_positions.extend(range(i1, i2))  # substituted or deleted reference bases
		edit_positions.extend(range(j1, j2))       # substituted or inserted edit bases
	return tuple(edit_positions), tuple(reference_positions)

def nucleotide_change_columns(reference_sequence, edit_sequence):
	"""Where the final allele differs from the reference, as two ';' separated index lists.

	The final allele carries the intended edit and, when one was made, the silent mutation,
	so this covers both. Two columns are needed because a deleted base has no position in
	the edited allele and an inserted base has none in the reference; for a substitution the
	two lists are identical. Indexes are 0-based into Reference_sequence and Edit_sequence.
	They are ';' separated because the output is written as naive comma-joined CSV.
	"""
	edit_positions, reference_positions = _nucleotide_change_indexes(str(reference_sequence), str(edit_sequence))
	return [';'.join(map(str, edit_positions)), ';'.join(map(str, reference_positions))]

##### Extract reference and edited sequence information
def process_sequence(input_sequence):

	# Check formatting is correct
	format_check = ''
	for i in input_sequence:
		if i == '(':
			format_check += '('
		elif i == ')':
			format_check += ')'
		elif i == '/':
			format_check += '/'
		elif i == '+':
			format_check += '+'
		elif i == '-':
			format_check += '-'

	# Check composition of input sequence
	if len(input_sequence) != sum([1 if x in ['A','T','C','G','(',')','+','-','/'] else 0 for x in input_sequence.upper()]):
		logger.error('Input sequence %s contains a character not in the following list: A,T,C,G,(,),+,-,/ ...' % str(input_sequence))
		sys.exit(1)

	# Check formatting
	if format_check.count('(') == format_check.count(')') and format_check.count('(') > 0: # Left and right parantheses equal
		if '((' not in format_check: # Checks both directions for nested parantheses
			if '()' not in format_check: # Checks for empty annotations
				if sum([1 if x in format_check else 0 for x in ['++','--','//','+-','+/','-+','-/','/+','/-','/(','+(','-(',')/',')+',')-']]) == 0:
					pass
				else:
					logger.error('Input sequence %s has more than one edit annotation per parantheses set (i.e. //,  +- , -/, etc.) ...' % str(input_sequence))
					sys.exit(1)
			else:
				logger.error('Input sequence %s has empty parantheses without an edit annotation (i.e. /,  + , -) ...' % str(input_sequence))
				sys.exit(1)
		else:
			logger.error('Input sequence %s has nested parantheses which is not allowed ...' % str(input_sequence))
			sys.exit(1)
	else:
		logger.error('Input sequence %s does not have full sets of parantheses ...' % str(input_sequence))
		sys.exit(1)

	# Create mapping between input format and reference and edit sequence
	editformat2sequence = {}
	edits = re.findall(r'\(.*?\)', input_sequence)
	for edit in edits:
		if '/' in edit:
			editformat2sequence[edit] = [edit.split('/')[0].replace('(',''), edit.split('/')[1].replace(')','')]
		elif '+' in edit:
			editformat2sequence[edit] = ['' , edit.split('+')[1].replace(')','')]
		elif '-' in edit:
			editformat2sequence[edit] = [edit.split('-')[1].replace(')',''), '']

	# Create mapping between edit number and reference and edit sequence
	editformat2sequence = {}
	editnumber2sequence = {}
	edit_idxs = [[m.start(), m.end()] for m in re.finditer(r'\(.*?\)', input_sequence)]
	edit_counter = 1
	edits_by_position = []  # (start, end, ref, alt, edit_counter), one entry per edit occurrence
	for edit_idx in edit_idxs:
		edit = input_sequence[edit_idx[0]:edit_idx[1]]

		# Create edit format and number to sequence map
		if '/' in edit:
			ref, alt = edit.split('/')[0].replace('(',''), edit.split('/')[1].replace(')','').lower()

		elif '+' in edit:
			ref, alt = '', edit.split('+')[1].replace(')','').lower()

		elif '-' in edit:
			ref, alt = edit.split('-')[1].replace(')',''), ''

		# Key on the edit text AND its start position so that two edit regions sharing identical
		# annotation text (e.g. two separate (A/T) substitutions at different positions) do not
		# overwrite each other in the dict.
		editformat2sequence['%s__at_%s' % (edit, str(edit_idx[0]))] = [ref, alt, edit_counter]
		editnumber2sequence[edit_counter] = [ref, alt]
		edits_by_position.append((edit_idx[0], edit_idx[1], ref, alt, edit_counter))

		edit_counter += 1

	edit_start = min([i.start() for i in re.finditer(r'\(', input_sequence)])
	edit_stop = max([i.start() for i in re.finditer(r'\)', input_sequence)])

	def _substitute_by_position(base_sequence, offset, value_index):
		# Rebuild base_sequence by substituting each parenthesized edit region (identified by its
		# absolute character span in input_sequence, shifted by `offset` for base_sequence) with
		# its ref (0), alt (1), or edit_counter (2) value. Edits are applied right-to-left so that
		# replacing one edit never shifts the positions of the edits still to be processed. This
		# replaces the previous str.replace(edit_text, ...) approach, which silently mis-numbered
		# (or dropped) edits that share identical annotation text at different positions.
		result = base_sequence
		for start, end, ref, alt, counter in sorted(edits_by_position, key = lambda e: e[0], reverse = True):
			value = (ref, alt, str(counter))[value_index]
			s, e = start - offset, end - offset
			result = result[:s] + value + result[e:]
		return(result)

	edit_span_sequence_w_ref = _substitute_by_position(input_sequence[edit_start:edit_stop + 1], edit_start, 0)
	edit_span_sequence_w_edit = _substitute_by_position(input_sequence[edit_start:edit_stop + 1], edit_start, 1)

	edit_start_in_ref = re.search(r'\(', input_sequence).start()
	edit_stop_in_ref_rev = re.search(r'\)', input_sequence[::-1]).start()

	edit_span_length_w_ref = len(edit_span_sequence_w_ref)
	edit_span_length_w_edit = len(edit_span_sequence_w_edit)

	reference_sequence = _substitute_by_position(input_sequence, 0, 0)
	edit_sequence = _substitute_by_position(input_sequence, 0, 1)
	editnumber_sequence = _substitute_by_position(input_sequence, 0, 2)

	return(editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev)

# Process sequence for saturating mutagenesis
def saturating_mutagenesis_input_sequences(target_name, target_sequence, sm_type):

	# Check formatting is correct
	format_check = ''
	for i in target_sequence:
		if i == '(':
			format_check += '('
		elif i == ')':
			format_check += ')'
		elif i == '/':
			format_check += '/'
		elif i == '+':
			format_check += '+'
		elif i == '-':
			format_check += '-'

	# Check for correct formatting of saturating mutagenesis input
	if len(target_sequence) != sum([1 if x in ['A','T','C','G', '(',')'] else 0 for x in target_sequence.upper()]):
		logger.error('Input sequence %s contains a character not in the following list: A,T,C,G,(,) ...' % str(target_sequence))
		sys.exit(1)

	# Check formatting
	if format_check.count('(') == format_check.count(')') and format_check.count('(') > 0: # Left and right parantheses equal
		if format_check.count('(') == 1:
			pass
		else:
			logger.error('Input sequence %s has more than one set of parantheses ...' % str(target_sequence))
			sys.exit(1)
	else:
		logger.error('Input sequence %s does not have full sets of parantheses ...' % str(target_sequence))
		sys.exit(1)

	parantheses_start = target_sequence.find('(')
	parantheses_stop = target_sequence.find(')')

	sequence_left = target_sequence[:parantheses_start]
	sequence_right = target_sequence[parantheses_stop + 1:]
	sequence_to_edit = target_sequence[parantheses_start + 1:parantheses_stop]

	sm_target_sequence_list = []
	sm_target_name_list = []
	sm_output_target_name_list = []
	sm_allowed_range_list = []
	sm_allowed_range = (len(sequence_left), len(sequence_left) + len(sequence_to_edit))

	def append_sm_target(output_target_name, sm_target_sequence, internal_suffix = None):
		"""Append a saturation-mutagenesis target.

		For non-base saturation mutagenesis modes, keep an internal target
		name that is unique so target_design entries do not overwrite each
		other. The output target name is recovered later by splitting at
		__internal_.
		"""
		if sm_type == 'base':
			internal_target_name = output_target_name
		else:
			internal_index = len(sm_target_name_list) + 1
			if internal_suffix is None:
				internal_suffix = str(internal_index)
			internal_target_name = '%s__internal_%s_%s' % (output_target_name, str(internal_index), internal_suffix)

		sm_target_name_list.append(internal_target_name)
		sm_target_sequence_list.append(sm_target_sequence)
		sm_output_target_name_list.append(output_target_name)
		sm_allowed_range_list.append(sm_allowed_range)

	def sorted_codons_by_difference(codon_ref, codon_edit_entry_list):
		"""Return all codons ranked by max nt difference, then codon usage."""
		return sorted(codon_edit_entry_list, key = lambda x: (-codon_nt_difference(codon_ref, x[0]), -x[1], x[0]))

	def sorted_codons_by_usage(codon_edit_entry_list):
		"""Return all codons ranked by codon usage."""
		return sorted(codon_edit_entry_list, key = lambda x: (-x[1], x[0]))

	if (sm_type == 'aa') | (sm_type == 'aa_subs'): # All Substitutions

		for base_index in range(0, len(sequence_to_edit), 3):

			codon_ref = sequence_to_edit[base_index:base_index + 3]

			if len(codon_ref) == 3:

				aa_ref = codon_dict[codon_ref][1]
				inner_sequence_left = sequence_to_edit[:base_index]
				inner_sequence_right = sequence_to_edit[base_index + 3:]

				aa_edit_list = [x for x in aa2codon if x != aa_ref]
				for aa_edit in aa_edit_list:

					codon_edit_list = sorted_codons_by_difference(codon_ref, aa2codon[aa_edit])
					for codon_edit_entry in codon_edit_list:

						codon_edit = codon_edit_entry[0]
						output_target_name = '%s_%s_%sto%s' % (target_name, str(int(base_index/3 + 1)), aa_ref, aa_edit)
						append_sm_target(output_target_name, sequence_left + inner_sequence_left + '(%s/%s)' % (codon_ref, codon_edit) + inner_sequence_right + sequence_right, '%s_%s' % (aa_edit, codon_edit))
	
	if (sm_type == 'aa') | (sm_type == 'aa_ins'): # All +1 Insertions
		
		for base_index in range(0, len(sequence_to_edit), 3):

			codon_ref = sequence_to_edit[base_index:base_index + 3]

			if len(codon_ref) == 3:

				aa_ref = codon_dict[codon_ref][1]
				inner_sequence_left = sequence_to_edit[:base_index]
				inner_sequence_right = sequence_to_edit[base_index+3:]

				aa_edit_list = [x for x in aa2codon if x != 'X']
				for aa_edit in aa_edit_list:

					codon_edit_list = sorted_codons_by_usage(aa2codon[aa_edit])
					for codon_edit_entry in codon_edit_list:

						codon_edit = codon_edit_entry[0]
						output_target_name = '%s_%s_%sto%s%s' % (target_name, str(int(base_index/3 + 1)), aa_ref, aa_ref, aa_edit)
						append_sm_target(output_target_name, sequence_left + inner_sequence_left + codon_ref + '(+%s)' % (codon_edit) + inner_sequence_right + sequence_right, '%s_%s' % (aa_edit, codon_edit))
	
	if (sm_type == 'aa') | (sm_type == 'aa_dels'): # All -1 Deletions
		
		for base_index in range(0, len(sequence_to_edit), 3):

			codon_ref = sequence_to_edit[base_index:base_index + 3]
			sequence_to_edit_and_right = sequence_to_edit + sequence_right
			codon2_ref = sequence_to_edit_and_right[base_index + 3:base_index + 6]

			if len(codon_ref) == 3:

				aa_ref = codon_dict[codon_ref][1]
				aa2_ref = codon_dict[codon2_ref][1]
				inner_sequence_left = sequence_to_edit[:base_index]
				inner_sequence_right = sequence_to_edit[base_index + 3:]

				output_target_name = '%s_%s_%s%sto%s' % (target_name, str(int(base_index/3 + 1)), aa_ref, aa2_ref, aa2_ref)
				append_sm_target(output_target_name, sequence_left + inner_sequence_left + '(-%s)' % (codon_ref) + inner_sequence_right + sequence_right, codon_ref)

	if sm_type == 'aa_silent': # All Silent Codon Substitutions

		for base_index in range(0, len(sequence_to_edit), 3):

			codon_ref = sequence_to_edit[base_index:base_index + 3]

			if len(codon_ref) == 3:

				aa_ref = codon_dict[codon_ref][1]
				inner_sequence_left = sequence_to_edit[:base_index]
				inner_sequence_right = sequence_to_edit[base_index + 3:]

				codon_edit_list = [x for x in aa2codon[aa_ref] if x[0] != codon_ref]
				codon_edit_list = sorted_codons_by_difference(codon_ref, codon_edit_list)
				for codon_edit_entry in codon_edit_list:

					codon_edit = codon_edit_entry[0]
					output_target_name = '%s_%s_%sto%s' % (target_name, str(int(base_index/3 + 1)), aa_ref, aa_ref)
					append_sm_target(output_target_name, sequence_left + inner_sequence_left + '(%s/%s)' % (codon_ref, codon_edit) + inner_sequence_right + sequence_right, codon_edit)

	if sm_type == 'base': # All Base Substitutions

		base_list = ['A','T','C','G']
		for base_index in range(len(sequence_to_edit)):

			base_ref = sequence_to_edit[base_index]
			inner_sequence_left = sequence_to_edit[:base_index]
			inner_sequence_right = sequence_to_edit[base_index + 1:]

			base_edit_list = [x for x in base_list if x != base_ref.upper()]
			for base_edit in base_edit_list:

				output_target_name = '%s_%s_%sto%s' % (target_name, str(base_index+1), base_ref, base_edit)
				append_sm_target(output_target_name, sequence_left + inner_sequence_left + '(%s/%s)' % (base_ref, base_edit) + inner_sequence_right + sequence_right)

	return(sm_target_name_list, sm_target_sequence_list, sm_output_target_name_list, sm_allowed_range_list)



def main():
	_configure_cli()

	##### Dictionary for to organize different DNA targets
	target_design = {}
	with open(file_in, 'r') as f:
		next(f) # Skip header line
	
		for line1 in f:
			line1 = line1.strip()

			if not line1: # Skip empty lines
				continue
		
			# Parse .txt files with space delimiter
			if file_in.lower().endswith('.txt'):
				parts = line1.split()

			# Parse .csv files with comma delimiter
			elif file_in.lower().endswith('.csv'):
				parts = line1.split(',')
		
			else:
				logger.error('Input file %s does not end with .txt or .csv ...' % str(file_in))
				sys.exit(1)
		
			# Accept either 2 or 3 columns
			if len(parts) < 2:
				logger.error(f"Line '{line1}' in {file_in} does not have at least 2 columns. ")
				sys.exit(1)

			target_name = parts[0].strip()
			target_sequence = parts[1].strip()

			# Optional third column: the amino acid number of the first residue the index
			# refers to. Under saturation mutagenesis that is the first residue of the
			# parenthesized region; otherwise it is the first residue of the target sequence.
			aa_index = 1
			if len(parts) > 2 and parts[2].strip():
				try:
					aa_index = int(parts[2].strip())
				except ValueError:
					logger.error("Third column '%s' in %s is not an integer amino acid index ..." % (parts[2].strip(), str(file_in)))
					sys.exit(1)

			target_sequence = target_sequence.upper()
		
			if saturation_mutagenesis:

				sm_target_name_list, sm_target_sequence_list, sm_output_target_name_list, sm_allowed_range_list = saturating_mutagenesis_input_sequences(target_name, target_sequence, saturation_mutagenesis)
			
				for sm_target_name, sm_target_sequence, sm_output_target_name, sm_allowed_range in zip(sm_target_name_list, sm_target_sequence_list, sm_output_target_name_list, sm_allowed_range_list):

					editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev = process_sequence(sm_target_sequence)

					# Initialize dictionary for the design of pegRNA spacers for each target sequence and intended edit(s)
					target_design[sm_target_name] = {'target_sequence':sm_target_sequence, 'editformat2sequence': editformat2sequence, 'editnumber2sequence': editnumber2sequence, 'reference_sequence': reference_sequence, 'edit_sequence': edit_sequence, 'editnumber_sequence': editnumber_sequence, 'edit_span_length': [edit_span_length_w_ref, edit_span_length_w_edit], 'edit_start_in_ref': edit_start_in_ref, 'edit_stop_in_ref_rev': edit_stop_in_ref_rev, 'pegRNA':{'+':[], '-':[]}, 'ngRNA':{'+':[], '-':[]}, 'output_target_name': sm_output_target_name, 'internal_target_name': sm_target_name, 'silent_mutation_allowed_range': sm_allowed_range, 'aa_index': aa_index}

			else:
				editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev = process_sequence(target_sequence)

				# Initialize dictionary for the design of pegRNA spacers for each target sequence and intended edit(s)
				target_design[target_name] = {'target_sequence':target_sequence, 'editformat2sequence': editformat2sequence, 'editnumber2sequence': editnumber2sequence, 'reference_sequence': reference_sequence, 'edit_sequence': edit_sequence, 'editnumber_sequence': editnumber_sequence, 'edit_span_length': [edit_span_length_w_ref, edit_span_length_w_edit], 'edit_start_in_ref': edit_start_in_ref, 'edit_stop_in_ref_rev': edit_stop_in_ref_rev, 'pegRNA':{'+':[], '-':[]}, 'ngRNA':{'+':[], '-':[]}, 'output_target_name': target_name, 'internal_target_name': target_name, 'silent_mutation_allowed_range': None, 'aa_index': aa_index}

	if len(target_design) == 0:
		logger.error('Input file %s does not have any entries. Make sure a column header is included (target_name,target_sequence) ...' % str(file_in))
		sys.exit(1)

	##### Find cut index and reformat PE format parameter
	if (pe_format.count('[') + pe_format.count(']')) == 2:

		if pe_format.count('/') == 1:

			# Find indices but shift when removing annotations
			cut_idx = re.search('/', pe_format).start()
			pam_start_idx = re.search(r'\[', pe_format).start()
			pam_end_idx = re.search(r'\]', pe_format).start()

			# Find pam and total PE format search length
			pam_length = pam_end_idx - pam_start_idx - 1
			pe_format_length = len(pe_format) - 3

			# Check if cut site is left of PAM
			if cut_idx < pam_start_idx:

				# Shift indices with removal of annotations
				pam_start_idx = pam_start_idx - 1
				pam_end_idx = pam_end_idx - 2
				spacer_start_idx = 0
				spacer_end_idx = pam_start_idx

			else:
				pam_end_idx = pam_end_idx - 1
				cut_idx = cut_idx - 2
				spacer_start_idx = pam_end_idx
				spacer_end_idx = len(pe_format) - 3
	
		else:
			logger.error('PE format parameter %s needs to cut site / within the spacer (i.e. NNNNNNNNNNNNNNNNN/NNN[NGG]) ...' % str(pe_format))
			sys.exit(1)

	else:
		logger.error('PE format parameter %s needs to have one [PAM] present in its sequence (i.e. NNNNNNNNNNNNNNNNN/NNN[NGG]) ...' % str(pe_format))
		sys.exit(1)

	# Remove annotations and convert into regex
	pe_format_rm_annotation = pe_format.replace('/', '').replace('[', '').replace(']', '')
	# print('---------- Prime editing spacer search parameters ----------')
	# print('PE format:\t%s' % pe_format_rm_annotation)
	# print('Spacer:\t\t%s' % pe_format_rm_annotation[spacer_start_idx:spacer_end_idx])
	# print('PAM:\t\t%s' % pe_format_rm_annotation[pam_start_idx:pam_end_idx])

	# Create PE format and PAM search sequences
	pe_format_search_plus = ''
	for base in pe_format_rm_annotation:
		pe_format_search_plus += iupac2bases(base)
	pe_format_search_minus = reverse_complement(pe_format_search_plus)

	pam_search = ''
	pam_sequence = pe_format_rm_annotation[pam_start_idx:pam_end_idx]
	for base in pam_sequence:
		pam_search += iupac2bases(base)

	# print('PE search (+):\t%s' % pe_format_search_plus)
	# print('PE search (-):\t%s' % pe_format_search_minus)
	# print('\n')

	##### Initialize data storage for output
	pe_design = {}
	logger.info('Searching for pegRNAs and nicking gRNAs for target sequences ...')
	counter = 1
	total_regions = len(target_design.keys())

	for target_name in target_design:

		# Store edit type
		edit_type = ''
		if "/" in target_design[target_name]['target_sequence']:
			edit_type += '& substitution'
		if "+" in target_design[target_name]['target_sequence']:
			edit_type += '& insertion'
		if "-" in target_design[target_name]['target_sequence']:
			edit_type += '& deletion'
		edit_type = edit_type[2:]

		# pegRNA spacer search for (+) and (-) strands with reference sequence
		reference_sequence = target_design[target_name]['reference_sequence']
		find_guides_ref_plus = [[m.start()] for m in re.finditer('(?=%s)' % pe_format_search_plus, reference_sequence, re.IGNORECASE)]
		find_guides_ref_minus = [[m.start()] for m in re.finditer('(?=%s)' % pe_format_search_minus, reference_sequence, re.IGNORECASE)]

		# pegRNA spacer search for (+) and (-) strands with edit number sequence
		editnumber_sequence = target_design[target_name]['editnumber_sequence']
		find_guides_editnumber_plus = [[m.start()] for m in re.finditer('(?=%s)' % pam_search.replace('[', '[123456789'), editnumber_sequence, re.IGNORECASE)]
		find_guides_editnumber_minus = [[m.start()] for m in re.finditer('(?=%s)' % reverse_complement(pam_search).replace('[', '[123456789'), editnumber_sequence, re.IGNORECASE)]

		editnumber2sequence = target_design[target_name]['editnumber2sequence']
		edit_sequence = target_design[target_name]['edit_sequence']

		# Find pegRNA spacers targeting (+) strand
		if find_guides_ref_plus:

			for match in find_guides_ref_plus:

				# Extract matched sequences and annotate type of prime editing
				full_search = reference_sequence[match[0]:match[0] + pe_format_length]
				spacer_sequence = full_search[spacer_start_idx:spacer_end_idx]
				extension_core_sequence = full_search[:cut_idx]
				downstream_sequence_ref = full_search[cut_idx:]
				downstream_sequence_length = len(downstream_sequence_ref)
				pam_ref = full_search[pam_start_idx:pam_end_idx]

				# Check to see if the extended non target strand is conserved in the edited strand
				try:
					extension_core_start_idx, extension_core_end_idx = re.search(extension_core_sequence, edit_sequence).start(), re.search(extension_core_sequence, edit_sequence).end()
					downstream_sequence_edit = edit_sequence[extension_core_end_idx:extension_core_end_idx + downstream_sequence_length]
					pam_edit = edit_sequence[extension_core_start_idx:extension_core_start_idx + pe_format_length][pam_start_idx:pam_end_idx]
				
					## Annotate pegRNA
					# Check if PAM is mutated relative to reference sequence
					if pam_ref == pam_edit.upper():
						pe_annotate = 'PAM_intact'

					else:
						# Check to see if mutation disrupts degenerate base positions within PAM
						if re.search(pam_search, pam_edit.upper()):
							pe_annotate = 'PAM_intact'

						else:
							pe_annotate = 'PAM_disrupted'

					# Store pegRNA spacer
					nick_ref_idx = match[0] + cut_idx
					nick_edit_idx = extension_core_start_idx + cut_idx
					target_design[target_name]['pegRNA']['+'].append([nick_ref_idx, nick_edit_idx, full_search, spacer_sequence, pam_ref, pam_edit, pe_annotate])

				except:
					continue

		# Find pegRNA spacers targeting (-) strand
		if find_guides_ref_minus:

			for match in find_guides_ref_minus:

				# Extract matched sequences and annotate type of prime editing
				full_search = reference_sequence[match[0]:match[0] + pe_format_length]
				spacer_sequence = full_search[pe_format_length - spacer_end_idx:pe_format_length - spacer_start_idx]
				extension_core_sequence = full_search[pe_format_length - cut_idx:]
				downstream_sequence_ref = full_search[:pe_format_length - cut_idx]
				downstream_sequence_length = len(downstream_sequence_ref)
				pam_ref = full_search[pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]

				# Check to see if the extended non target strand is conserved in the edited strand
				try:
					extension_core_start_idx, extension_core_end_idx = re.search(extension_core_sequence, edit_sequence).start(), re.search(extension_core_sequence, edit_sequence).end()
					downstream_sequence_edit = edit_sequence[extension_core_start_idx - downstream_sequence_length:extension_core_start_idx]
					pam_edit = edit_sequence[extension_core_end_idx - pe_format_length:extension_core_end_idx][pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]
				
					## Annotate pegRNA
					# Check if PAM is mutated relative to reference sequence
					if pam_ref == pam_edit.upper():
						pe_annotate = 'PAM_intact'

					else:
						# Check to see if mutation disrupts degenerate base positions within PAM
						if re.search(reverse_complement(pam_search), pam_edit.upper()):
							pe_annotate = 'PAM_intact'

						else:
							pe_annotate = 'PAM_disrupted'

					# Store pegRNA spacer
					nick_ref_idx = match[0] + (pe_format_length - cut_idx)
					nick_edit_idx = extension_core_start_idx - downstream_sequence_length + (pe_format_length - cut_idx)
					target_design[target_name]['pegRNA']['-'].append([nick_ref_idx, nick_edit_idx, full_search, spacer_sequence, pam_ref, pam_edit, pe_annotate])

				except:
					continue

		# Find ngRNA spacers targeting (+) strand
		if find_guides_editnumber_plus:

			for match in find_guides_editnumber_plus:

				# Extract matched sequences and annotate type of prime editing
				full_search = editnumber_sequence[:match[0] + pam_length]
			
				full_search2ref = full_search
				full_search2edit = full_search
				for edit_number in editnumber2sequence:
					full_search2ref = full_search2ref.replace(str(edit_number), editnumber2sequence[edit_number][0])
					full_search2edit = full_search2edit.replace(str(edit_number), editnumber2sequence[edit_number][1])

				if len(full_search2edit[-pe_format_length:]) == pe_format_length:

					# Identify ngRNA sequence information from edit sequence
					full_search_edit = full_search2edit[-pe_format_length:]
					spacer_sequence_edit = full_search_edit[spacer_start_idx:spacer_end_idx]
					pam_edit = full_search_edit[pam_start_idx:pam_end_idx]

					# Use reference sequence to find nick index
					full_search_ref = full_search2ref[-pe_format_length:]
					spacer_sequence_ref = full_search_ref[spacer_start_idx:spacer_end_idx]
					pam_ref = full_search_ref[pam_start_idx:pam_end_idx]

					# Annotate ngRNA
					if spacer_sequence_edit.upper()	== spacer_sequence_ref.upper():
						ng_annotate = 'PE3'
					else:
						if spacer_sequence_edit.upper()[-10:] == spacer_sequence_ref.upper()[-10:]:
							ng_annotate = 'PE3b-nonseed'
						else:
							ng_annotate = 'PE3b-seed'

					# Store ngRNA spacer
					nick_ref_idx = re.search(full_search_ref, reference_sequence).end() - (pe_format_length - cut_idx)
					nick_edit_start_idx = re.search(spacer_sequence_edit, edit_sequence).start()
					nick_edit_end_idx = re.search(spacer_sequence_edit, edit_sequence).end()
					target_design[target_name]['ngRNA']['+'].append([nick_ref_idx, nick_edit_start_idx, nick_edit_end_idx, full_search_edit, spacer_sequence_edit, pam_edit, ng_annotate])

		# Find ngRNA spacers targeting (-) strand
		if find_guides_editnumber_minus:

			for match in find_guides_editnumber_minus:

				# Extract matched sequences and annotate type of prime editing
				full_search = editnumber_sequence[match[0]:]
			
				full_search2ref = full_search
				full_search2edit = full_search
				for edit_number in editnumber2sequence:
					full_search2ref = full_search2ref.replace(str(edit_number), editnumber2sequence[edit_number][0])
					full_search2edit = full_search2edit.replace(str(edit_number), editnumber2sequence[edit_number][1])

				if len(full_search2edit[:pe_format_length]) == pe_format_length:

					# Identify ngRNA sequence information from edit sequence
					full_search_edit = full_search2edit[:pe_format_length]
					spacer_sequence_edit = full_search_edit[pe_format_length - spacer_end_idx:pe_format_length - spacer_start_idx]
					pam_edit = full_search_edit[pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]

					# Use reference sequence to find nick index
					full_search_ref = full_search2ref[:pe_format_length]
					spacer_sequence_ref = full_search_ref[pe_format_length - spacer_end_idx:pe_format_length - spacer_start_idx]
					pam_ref = full_search_ref[pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]

					# Annotate ngRNA
					if spacer_sequence_edit.upper()	== spacer_sequence_ref.upper():
						ng_annotate = 'PE3'
					else:
						if spacer_sequence_edit.upper()[:10] == spacer_sequence_ref.upper()[:10]:
							ng_annotate = 'PE3b-nonseed'
						else:
							ng_annotate = 'PE3b-seed'

					# Store ngRNA spacer
					nick_ref_idx = re.search(full_search_ref, reference_sequence).start() + (pe_format_length - cut_idx)
					nick_edit_start_idx = re.search(spacer_sequence_edit, edit_sequence).start()
					nick_edit_end_idx = re.search(spacer_sequence_edit, edit_sequence).end()
					target_design[target_name]['ngRNA']['-'].append([nick_ref_idx, nick_edit_start_idx, nick_edit_end_idx, full_search_edit, spacer_sequence_edit, pam_edit, ng_annotate])

		# Grab index information of edits to introduce to target sequence
		edit_start_in_ref = int(target_design[target_name]['edit_start_in_ref'])
		edit_stop_in_ref_rev = int(target_design[target_name]['edit_stop_in_ref_rev'])
		edit_span_length_w_ref = int(target_design[target_name]['edit_span_length'][0])
		edit_span_length_w_edit = int(target_design[target_name]['edit_span_length'][1])

		# Initialize pegRNA and ngRNA design dictionary
		pe_design[target_name] = {}

		### Separate genome-wide, saturation mutagenesis, and regular pegRNA design here

		# Design for genome-wide or saturation mutagenesis screening applications
		if genome_wide_design or saturation_mutagenesis:
		
			# Design pegRNAs targeting the (+) strand
			for peg_plus in target_design[target_name]['pegRNA']['+']:

				pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_plus

				pe_annotate_constant = pe_annotate

				# See if pegRNA spacer can introduce all edits and correct orientation (downstream)
				nick2edit_length = edit_start_in_ref - pe_nick_ref_idx
				if nick2edit_length >= 0:

					# See if RTT length can reach entire edit with homology downstream constraint
					silent_mutation_edit_sequence = ''
					nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
					rtt_length = nick2lastedit_length + homology_downstream
					if rtt_length < rtt_max_length_pooled:

						pbs_length = pbs_length_pooled
						pe_pam_ref_silent_mutation = ''

						# Construct pegRNA extension to encode intended edit(s)
						# Silent mutations only work for NGG PAMs
						pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])
						pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_max_length_pooled])
						nick_aa_index = int(pe_nick_edit_idx)%3
						silent_mutation_relative_to_edit = NO_SILENT_MUTATION
						silent_mutation_selected = None
						silent_options = [None]
						if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]'): # Generate silent mutations

							# The intended edit may itself have altered the PAM. ranked_silent_mutations() reads the
							# PAM out of the edited allele, so a silent codon only earns PAM-disruption credit when the
							# edit left the PAM intact; when the edit already broke it, the search falls through to the
							# next criterion and still adds a separate silent mutation inside the allowed region.
							edit_start_abs = pe_nick_edit_idx + nick2edit_length
							edit_stop_abs = edit_start_abs + edit_span_length_w_edit
							# Make room if the arm leaves none: a longer RTT pushes the arm 3'
							rtt_length, silent_candidates = extend_rtt_for_silent_mutations(
								target_design[target_name], edit_sequence, pe_nick_edit_idx, rtt_length,
								'+', pam_search, edit_start_abs, edit_stop_abs, homology_downstream,
								rtt_max_length_pooled)
							pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

							# Saturation keeps its single best silent mutation; the no-silent design is offered too
							silent_options = with_no_silent_mutation_option(silent_candidates[:1])

						base_pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])
						base_pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_max_length_pooled])
						# The ngRNA block below reads pegid and silent_mutation_edit_sequence, so the first
						# option is restored once the loop ends to leave its behavior unchanged.
						primary_pegid = None
						primary_silent_mutation_edit_sequence = ''
						for silent_mutation_selected in silent_options:

							pe_pam_ref_silent_mutation = ''
							silent_mutation_relative_to_edit = NO_SILENT_MUTATION
							silent_mutation_edit_sequence = ''
							pegRNA_ext = base_pegRNA_ext
							pegRNA_ext_max = base_pegRNA_ext_max

							if silent_mutation_selected is not None:
								silent_mutation_edit_sequence = silent_mutation_selected['silent_mutation_edit_sequence']
								pegRNA_ext = reverse_complement(silent_mutation_edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])
								pegRNA_ext_max = reverse_complement(silent_mutation_edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_max_length_pooled])
								pe_silent_mutation = silent_mutation_selected['codon_ref'] + '-to-' + silent_mutation_selected['codon_edit'].lower()

								if silent_mutation_selected['disrupts_pam']:
									pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + silent_mutation_edit_sequence[pe_nick_edit_idx + 3:pe_nick_edit_idx + 6].lower()


								if silent_mutation_selected['codon_start_abs'] >= edit_stop_abs:
									silent_mutation_relative_to_edit = 'downstream'
								else:
									silent_mutation_relative_to_edit = 'upstream'
								if not primary_silent_mutation_edit_sequence:
									primary_silent_mutation_edit_sequence = silent_mutation_edit_sequence

							# Classify the final allele: what can this pegRNA no longer re-cut?
							pe_annotate = pam_spacer_annotation(reference_sequence,
								silent_mutation_edit_sequence if silent_mutation_edit_sequence else edit_sequence,
								pe_nick_ref_idx, pe_nick_edit_idx, '+', pam_search,
								cut_idx, spacer_end_idx, pam_start_idx, pam_end_idx)
							# Enforce saturation-mutagenesis silent mutation boundaries.
							# If the extra silent mutation falls outside the parenthesized region,
							# revert to the original edit sequence and original PAM annotation.
							if saturation_mutagenesis and len(silent_mutation_edit_sequence) > 0 and not silent_mutation_edit_sequence_in_allowed_range(target_design[target_name], edit_sequence, silent_mutation_edit_sequence):
								pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])
								pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_max_length_pooled])
								pe_pam_ref_silent_mutation = ''
								pe_annotate = pe_annotate_constant
								silent_mutation_edit_sequence = ''
								silent_mutation_relative_to_edit = NO_SILENT_MUTATION

							# Create pegRNA ID based on annotation
							if pe_annotate == 'PAM_disrupted':
								pe_annotate_code = 0
							elif pe_annotate == 'Spacer_disrupted':
								pe_annotate_code = 1
							else:
								pe_annotate_code = 2

							pegid = '_'.join(map(str, [str(pe_annotate_code) + '0'*(3 - len(str(abs(nick2lastedit_length)))) + str(abs(nick2lastedit_length)), pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+']))
					
							# Check to see if pegRNA extension is within input sequence
							# ngRNAs for this allele: a PE3b spacer is cut from the sequence this pegRNA
							# installs, so it is built per silent mutation option rather than shared
							ngRNA_entries = build_ngRNA_entries(target_design[target_name], '-',
								silent_mutation_edit_sequence, pe_nick_ref_idx, silent_mutation, pe_format,
								nicking_distance_minimum, nicking_distance_maximum, nicking_distance_pooled,
								ranked = True)

							if len(pegRNA_ext) == (pbs_length + rtt_length):

								# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
								if pegid not in pe_design[target_name]:

									# First list is for peg extension, second list is for nicking guide
									pe_design[target_name][pegid] = [[],[]]
						
								# Store pegRNA design
								if silent_mutation_edit_sequence == '':
									if pe_pam_ref_silent_mutation == '':
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])

									else:
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref_silent_mutation, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])
								else:
									if pe_pam_ref_silent_mutation == '':
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])

									else:
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref_silent_mutation, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])

							if primary_pegid is None:
								primary_pegid = pegid

						silent_mutation_edit_sequence = primary_silent_mutation_edit_sequence
						if primary_pegid is not None:
							pegid = primary_pegid


			# Design pegRNAs targeting the (-) strand
			for peg_minus in target_design[target_name]['pegRNA']['-']:

				pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_minus

				pe_annotate_constant = pe_annotate

				# See if pegRNA spacer can introduce all edits
				nick2edit_length = edit_stop_in_ref_rev - (len(reference_sequence) - pe_nick_ref_idx)
				if nick2edit_length >= 0:

					# See if RT length can reach entire edit
					silent_mutation_edit_sequence = ''
					nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
					rtt_length = nick2lastedit_length + homology_downstream
					if rtt_length < rtt_max_length_pooled:

						pbs_length = pbs_length_pooled
						pe_pam_ref_silent_mutation = ''

						# Construct pegRNA extension to encode intended edit(s)
						# Silent mutations only work for NGG PAMs
						pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
						pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + pbs_length]
						nick_aa_index = int(pe_nick_edit_idx)%3
						silent_mutation_relative_to_edit = NO_SILENT_MUTATION
						silent_mutation_selected = None
						silent_options = [None]
						if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]'):

							# Mirror of the (+) strand case: the RTT runs 5' of the nick, so the PAM sits at
							# [nick-6, nick-3) and the edit lies 5' of the nick. As on the (+) strand, a silent
							# codon earns PAM-disruption credit only when the intended edit left the PAM intact.
							edit_stop_abs = pe_nick_edit_idx - nick2edit_length
							edit_start_abs = pe_nick_edit_idx - nick2lastedit_length
							# Make room if the arm leaves none: a longer RTT pushes the arm 3'
							rtt_length, silent_candidates = extend_rtt_for_silent_mutations(
								target_design[target_name], edit_sequence, pe_nick_edit_idx, rtt_length,
								'-', pam_search, edit_start_abs, edit_stop_abs, homology_downstream,
								rtt_max_length_pooled)
							pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]

							# Saturation keeps its single best silent mutation; the no-silent design is offered too
							silent_options = with_no_silent_mutation_option(silent_candidates[:1])

						base_pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
						base_pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + pbs_length]
						# The ngRNA block below reads pegid and silent_mutation_edit_sequence, so the first
						# option is restored once the loop ends to leave its behavior unchanged.
						primary_pegid = None
						primary_silent_mutation_edit_sequence = ''
						for silent_mutation_selected in silent_options:

							pe_pam_ref_silent_mutation = ''
							silent_mutation_relative_to_edit = NO_SILENT_MUTATION
							silent_mutation_edit_sequence = ''
							pegRNA_ext = base_pegRNA_ext
							pegRNA_ext_max = base_pegRNA_ext_max

							if silent_mutation_selected is not None:
								silent_mutation_edit_sequence = silent_mutation_selected['silent_mutation_edit_sequence']
								pegRNA_ext = silent_mutation_edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
								pegRNA_ext_max = silent_mutation_edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + pbs_length]
								pe_silent_mutation = reverse_complement(silent_mutation_selected['codon_ref']) + '-to-' + reverse_complement(silent_mutation_selected['codon_edit']).lower()

								if silent_mutation_selected['disrupts_pam']:
									pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(silent_mutation_edit_sequence[pe_nick_edit_idx - 6:pe_nick_edit_idx - 3]).lower()


								if silent_mutation_selected['codon_start_abs'] >= edit_stop_abs:
									silent_mutation_relative_to_edit = 'downstream'
								else:
									silent_mutation_relative_to_edit = 'upstream'
								if not primary_silent_mutation_edit_sequence:
									primary_silent_mutation_edit_sequence = silent_mutation_edit_sequence

							# Classify the final allele: what can this pegRNA no longer re-cut?
							pe_annotate = pam_spacer_annotation(reference_sequence,
								silent_mutation_edit_sequence if silent_mutation_edit_sequence else edit_sequence,
								pe_nick_ref_idx, pe_nick_edit_idx, '-', pam_search,
								cut_idx, spacer_end_idx, pam_start_idx, pam_end_idx)
							# Enforce saturation-mutagenesis silent mutation boundaries.
							# If the extra silent mutation falls outside the parenthesized region,
							# revert to the original edit sequence and original PAM annotation.
							if saturation_mutagenesis and len(silent_mutation_edit_sequence) > 0 and not silent_mutation_edit_sequence_in_allowed_range(target_design[target_name], edit_sequence, silent_mutation_edit_sequence):
								pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
								pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + pbs_length]
								pe_pam_ref_silent_mutation = ''
								pe_annotate = pe_annotate_constant
								silent_mutation_edit_sequence = ''
								silent_mutation_relative_to_edit = NO_SILENT_MUTATION

							# Create pegRNA ID based on annotation
							if pe_annotate == 'PAM_disrupted':
								pe_annotate_code = 0
							elif pe_annotate == 'Spacer_disrupted':
								pe_annotate_code = 1
							else:
								pe_annotate_code = 2

							pegid = '_'.join(map(str, [str(pe_annotate_code) + '0'*(3 - len(str(abs(nick2lastedit_length)))) + str(abs(nick2lastedit_length)), pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-']))
					
							# Check to see if pegRNA extension is within input sequence
							# ngRNAs for this allele: a PE3b spacer is cut from the sequence this pegRNA
							# installs, so it is built per silent mutation option rather than shared
							ngRNA_entries = build_ngRNA_entries(target_design[target_name], '+',
								silent_mutation_edit_sequence, pe_nick_ref_idx, silent_mutation, pe_format,
								nicking_distance_minimum, nicking_distance_maximum, nicking_distance_pooled,
								ranked = True)

							if len(pegRNA_ext) == (pbs_length + rtt_length):

								# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
								if pegid not in pe_design[target_name]:

									# First list is for peg extension, second list is for nicking guide
									pe_design[target_name][pegid] = [[],[]]
						
								# Store pegRNA design
								if silent_mutation_edit_sequence == '':
									if pe_pam_ref_silent_mutation == '':
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), reverse_complement(pe_pam_ref), pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])
							
									else:
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), pe_pam_ref_silent_mutation, pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])
						
								else:
									if pe_pam_ref_silent_mutation == '':
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), reverse_complement(pe_pam_ref), pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])
							
									else:
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), pe_pam_ref_silent_mutation, pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])

							if primary_pegid is None:
								primary_pegid = pegid

						silent_mutation_edit_sequence = primary_silent_mutation_edit_sequence
						if primary_pegid is not None:
							pegid = primary_pegid


			# Sort pegRNAs and ngRNAs and filter for top designs
			pe_design[target_name] = dict(sorted(pe_design[target_name].items(), key=lambda v: int(v[0].split('_')[0])))

			if counter%1000 == 0:
				logger.info('Completed pegRNA and ngRNA search for %s out of %s sites ...' % (counter, total_regions))
			counter += 1

		# Normal pegRNA design with PBS and RTT parameter ranges
		else:
	
			# Design pegRNAs targeting the (+) strand
			for peg_plus in target_design[target_name]['pegRNA']['+']:

				pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_plus
				pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+']))

				pe_annotate_constant = pe_annotate

				# See if pegRNA spacer can introduce all edits
				nick2edit_length = edit_start_in_ref - pe_nick_ref_idx
				if nick2edit_length >= 0:

					# Loop through RTT lengths
					silent_mutation_edit_sequence = ''
					for rtt_length in rtt_length_list:

						# See if RT length can reach entire edit
						nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
						if nick2lastedit_length + homology_downstream < rtt_length:

							# Loop through PBS lengths
							for pbs_length in pbs_length_list:
								pe_pam_ref_silent_mutation = ''

								# Construct pegRNA extension to encode intended edit(s)
								# Silent mutations only work for NGG PAMs
								pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])
								nick_aa_index = int(pe_nick_edit_idx)%3
								silent_mutation_relative_to_edit = NO_SILENT_MUTATION
								silent_options = [None]
								if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]'):

									# Non-saturation (+) strand path. silent_mutation_allowed_range is None here, so every
									# codon inside the RTT is eligible; ranking is unchanged: PAM, then spacer, then
									# nucleotide difference, then codon usage.
									edit_start_abs = pe_nick_edit_idx + nick2edit_length
									edit_stop_abs = edit_start_abs + edit_span_length_w_edit
									silent_candidates = ranked_silent_mutations(target_design[target_name], edit_sequence,
										pe_nick_edit_idx, rtt_length, '+', pam_search, edit_start_abs, edit_stop_abs, homology_downstream)


									silent_options = with_no_silent_mutation_option(
										select_silent_mutations(silent_candidates, number_of_silent_mutations))

								# The ngRNA block after these loops reads silent_mutation_edit_sequence, so the best
								# option is restored once the loop ends to leave that behavior unchanged.
								primary_silent_mutation_edit_sequence = ''
								base_pegRNA_ext = pegRNA_ext
								for silent_mutation_selected in silent_options:

									pe_pam_ref_silent_mutation = ''
									silent_mutation_relative_to_edit = NO_SILENT_MUTATION
									silent_mutation_edit_sequence = ''
									pegRNA_ext = base_pegRNA_ext

									if silent_mutation_selected is not None:
										silent_mutation_edit_sequence = silent_mutation_selected['silent_mutation_edit_sequence']
										pegRNA_ext = reverse_complement(silent_mutation_edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])
										pe_silent_mutation = silent_mutation_selected['codon_ref'] + '-to-' + silent_mutation_selected['codon_edit'].lower()

										if silent_mutation_selected['disrupts_pam']:
											pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + silent_mutation_edit_sequence[pe_nick_edit_idx + 3:pe_nick_edit_idx + 6].lower()


										if silent_mutation_selected['codon_start_abs'] >= edit_stop_abs:
											silent_mutation_relative_to_edit = 'downstream'
										else:
											silent_mutation_relative_to_edit = 'upstream'
										if not primary_silent_mutation_edit_sequence:
											primary_silent_mutation_edit_sequence = silent_mutation_edit_sequence

									# Classify the final allele: what can this pegRNA no longer re-cut?
									pe_annotate = pam_spacer_annotation(reference_sequence,
										silent_mutation_edit_sequence if silent_mutation_edit_sequence else edit_sequence,
										pe_nick_ref_idx, pe_nick_edit_idx, '+', pam_search,
										cut_idx, spacer_end_idx, pam_start_idx, pam_end_idx)
									# Check to see if pegRNA extension is within input sequence
									# ngRNAs for this allele: a PE3b spacer is cut from the sequence this pegRNA
									# installs, so it is built per silent mutation option rather than shared
									ngRNA_entries = build_ngRNA_entries(target_design[target_name], '-',
										silent_mutation_edit_sequence, pe_nick_ref_idx, silent_mutation, pe_format,
										nicking_distance_minimum, nicking_distance_maximum, nicking_distance_pooled,
										ranked = False)

									if len(pegRNA_ext) == (pbs_length + rtt_length):

										# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
										if pegid not in pe_design[target_name]:

											# First list is for peg extension, second list is for nicking guide
											pe_design[target_name][pegid] = [[],[]]

										# Store pegRNA design
										if silent_mutation_edit_sequence == '':
											if pe_pam_ref_silent_mutation == '':
												pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])

											else:
												pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref_silent_mutation, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])
										else:
											if pe_pam_ref_silent_mutation == '':
												pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])

											else:
												pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref_silent_mutation, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])

								silent_mutation_edit_sequence = primary_silent_mutation_edit_sequence


			# Design pegRNAs targeting the (-) strand
			for peg_minus in target_design[target_name]['pegRNA']['-']:

				pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_minus
				pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-']))

				pe_annotate_constant = pe_annotate

				# See if pegRNA spacer can introduce all edits
				nick2edit_length = edit_stop_in_ref_rev - (len(reference_sequence) - pe_nick_ref_idx)
				if nick2edit_length >= 0:

					# Loop through RTT lengths
					silent_mutation_edit_sequence = ''
					for rtt_length in rtt_length_list:

						# See if RT length can reach entire edit
						nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
						if nick2lastedit_length + homology_downstream < rtt_length:

							# Loop through PBS lengths
							for pbs_length in pbs_length_list:
								pe_pam_ref_silent_mutation = ''

								# Construct pegRNA extension to encode intended edit(s)
								# Silent mutations only work for NGG PAMs
								pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
								nick_aa_index = int(pe_nick_edit_idx)%3
								silent_mutation_relative_to_edit = NO_SILENT_MUTATION
								silent_options = [None]
								if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]'):

									# Non-saturation (-) strand path: RTT runs 5' of the nick, PAM at [nick-6, nick-3).
									edit_stop_abs = pe_nick_edit_idx - nick2edit_length
									edit_start_abs = pe_nick_edit_idx - nick2lastedit_length
									silent_candidates = ranked_silent_mutations(target_design[target_name], edit_sequence,
										pe_nick_edit_idx, rtt_length, '-', pam_search, edit_start_abs, edit_stop_abs, homology_downstream)


									silent_options = with_no_silent_mutation_option(
										select_silent_mutations(silent_candidates, number_of_silent_mutations))

								# The ngRNA block after these loops reads silent_mutation_edit_sequence, so the best
								# option is restored once the loop ends to leave that behavior unchanged.
								primary_silent_mutation_edit_sequence = ''
								base_pegRNA_ext = pegRNA_ext
								for silent_mutation_selected in silent_options:

									pe_pam_ref_silent_mutation = ''
									silent_mutation_relative_to_edit = NO_SILENT_MUTATION
									silent_mutation_edit_sequence = ''
									pegRNA_ext = base_pegRNA_ext

									if silent_mutation_selected is not None:
										silent_mutation_edit_sequence = silent_mutation_selected['silent_mutation_edit_sequence']
										pegRNA_ext = silent_mutation_edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
										pe_silent_mutation = reverse_complement(silent_mutation_selected['codon_ref']) + '-to-' + reverse_complement(silent_mutation_selected['codon_edit']).lower()

										if silent_mutation_selected['disrupts_pam']:
											pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(silent_mutation_edit_sequence[pe_nick_edit_idx - 6:pe_nick_edit_idx - 3]).lower()


										if silent_mutation_selected['codon_start_abs'] >= edit_stop_abs:
											silent_mutation_relative_to_edit = 'downstream'
										else:
											silent_mutation_relative_to_edit = 'upstream'
										if not primary_silent_mutation_edit_sequence:
											primary_silent_mutation_edit_sequence = silent_mutation_edit_sequence

									# Classify the final allele: what can this pegRNA no longer re-cut?
									pe_annotate = pam_spacer_annotation(reference_sequence,
										silent_mutation_edit_sequence if silent_mutation_edit_sequence else edit_sequence,
										pe_nick_ref_idx, pe_nick_edit_idx, '-', pam_search,
										cut_idx, spacer_end_idx, pam_start_idx, pam_end_idx)
									# Check to see if pegRNA extension is within input sequence
									# ngRNAs for this allele: a PE3b spacer is cut from the sequence this pegRNA
									# installs, so it is built per silent mutation option rather than shared
									ngRNA_entries = build_ngRNA_entries(target_design[target_name], '+',
										silent_mutation_edit_sequence, pe_nick_ref_idx, silent_mutation, pe_format,
										nicking_distance_minimum, nicking_distance_maximum, nicking_distance_pooled,
										ranked = False)

									if len(pegRNA_ext) == (pbs_length + rtt_length):

										# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
										if pegid not in pe_design[target_name]:

											# First list is for peg extension, second list is for nicking guide
											pe_design[target_name][pegid] = [[],[]]

										# Store pegRNA design
										if silent_mutation_edit_sequence == '':
											if pe_pam_ref_silent_mutation == '':
												pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), reverse_complement(pe_pam_ref), pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])
									
											else:
												pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), pe_pam_ref_silent_mutation, pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])
								
										else:
											if pe_pam_ref_silent_mutation == '':
												pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), reverse_complement(pe_pam_ref), pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])
									
											else:
												pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), pe_pam_ref_silent_mutation, pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries])

								silent_mutation_edit_sequence = primary_silent_mutation_edit_sequence


			if counter%1000 == 0:
				logger.info('Completed pegRNA and ngRNA search for %s out of %s sites ...' % (counter, total_regions))
			counter += 1

	logger.info('Completed pegRNA and ngRNA search for %s out of %s sites ...' % (counter - 1, total_regions))

	# Output pegRNAs
	pegRNAs_summary_f = '%s_PrimeDesign.csv' % str(time.strftime("%Y%m%d_%I.%M.%S", time.localtime()))
	logger.info('Writing pegRNA and ngRNA designs into output file %s ...' % pegRNAs_summary_f)

	counter = 1
	with open(out_dir + '/%s' % pegRNAs_summary_f, 'w') as f:

		f.write(','.join(map(str, ['Target_name', 'Target_sequence', 'pegRNA_number', 'gRNA_type', 'Spacer_sequence', 'Spacer_GC_content', 'PAM_sequence', 'Extension_sequence', 'Strand', 'Annotation', 'pegRNA-to-edit_distance', 'Nick_index', 'ngRNA-to-pegRNA_distance', 'PBS_length', 'PBS_GC_content', 'RTT_length', 'RTT_GC_content', 'First_extension_nucleotide', 'Spacer_sequence_order_TOP', 'Spacer_sequence_order_BOTTOM', 'pegRNA_extension_sequence_order_TOP', 'pegRNA_extension_sequence_order_BOTTOM', 'Edit_type', 'Reference_sequence', 'Edit_sequence','Silent_mutation_relative_to_edit'] + SILENT_CODON_COLUMNS + NUCLEOTIDE_CHANGE_COLUMNS)) + '\n')
		for target_name in pe_design:

			if genome_wide_design or saturation_mutagenesis:

				if filter_c1_extension:
				
					peg_count = 0
					for pegid in list(pe_design[target_name].keys()):

						ng_continue = True

						# Write pegRNAs
						for pegRNA_entry in pe_design[target_name][pegid][0]:
							pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, pe_strand, pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries = pegRNA_entry

							pegRNA_ext_first_base = pegRNA_ext[0]
							spacer_gc_content = gc_content(pe_spacer_sequence)
							pbs_gc_content = gc_content(pegRNA_ext[rtt_length:])
							rtt_gc_content = gc_content(pegRNA_ext[:rtt_length])

							if pegRNA_ext_first_base.upper() == 'C':

								# Find minimum non-C index extending the RTT template
								shift_rtt_index_list = [pegRNA_ext_max[:-len(pegRNA_ext)][::-1].upper().find('A'), pegRNA_ext_max[:-len(pegRNA_ext)][::-1].upper().find('G'), pegRNA_ext_max[:-len(pegRNA_ext)][::-1].upper().find('T')]
								shift_rtt_index_list = [x for x in shift_rtt_index_list if x != -1]

								# Make sure there are non-C indices
								if len(shift_rtt_index_list) > 0:
									shift_rtt_index = min(shift_rtt_index_list) + 1
									rtt_length += shift_rtt_index

									pegRNA_ext = pegRNA_ext_max[-len(pegRNA_ext) - shift_rtt_index:]
									pegRNA_ext_first_base = pegRNA_ext[0]
									spacer_gc_content = gc_content(pe_spacer_sequence)
									pbs_gc_content = gc_content(pegRNA_ext[rtt_length:])
									rtt_gc_content = gc_content(pegRNA_ext[:rtt_length])

									if pe_spacer_sequence[0].upper() == 'G':
										spacer_oligo_top = 'cacc' + pe_spacer_sequence + 'gtttt'
										spacer_oligo_bottom = 'ctctaaaac' + reverse_complement(pe_spacer_sequence)

									else:
										spacer_oligo_top = 'caccG' + pe_spacer_sequence + 'gtttt'
										spacer_oligo_bottom = 'ctctaaaac' + reverse_complement('G' + pe_spacer_sequence)

									pegext_oligo_top = 'gtgc' + pegRNA_ext
									pegext_oligo_bottom = 'aaaa' + reverse_complement(pegRNA_ext)

									if filter_homopolymer_ts:

										if 'TTTT' not in pe_spacer_sequence:
											f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + codon_columns(target_design[target_name], silent_mutation_selected) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')
											peg_count += 1

										else:
											ng_continue = False

									else:

										f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + codon_columns(target_design[target_name], silent_mutation_selected) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')
										peg_count += 1

							else:

								if pe_spacer_sequence[0].upper() == 'G':
									spacer_oligo_top = 'cacc' + pe_spacer_sequence + 'gtttt'
									spacer_oligo_bottom = 'ctctaaaac' + reverse_complement(pe_spacer_sequence)

								else:
									spacer_oligo_top = 'caccG' + pe_spacer_sequence + 'gtttt'
									spacer_oligo_bottom = 'ctctaaaac' + reverse_complement('G' + pe_spacer_sequence)

								pegext_oligo_top = 'gtgc' + pegRNA_ext
								pegext_oligo_bottom = 'aaaa' + reverse_complement(pegRNA_ext)

								if filter_homopolymer_ts:

									if 'TTTT' not in pe_spacer_sequence:
										f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + codon_columns(target_design[target_name], silent_mutation_selected) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')
										peg_count += 1
									
									else:
										ng_continue = False

								else:

									f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + codon_columns(target_design[target_name], silent_mutation_selected) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')
									peg_count += 1

							# Write ngRNAs
							if ng_continue:
								for ngRNA_entry in ngRNA_entries[:number_of_ngrnas]:
									ng_code, ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, ng_strand, nick_distance = ngRNA_entry

									spacer_gc_content = gc_content(ng_spacer_sequence_edit)

									if ng_spacer_sequence_edit[0].upper() == 'G':
										spacer_oligo_top = 'cacc' + ng_spacer_sequence_edit
										spacer_oligo_bottom = 'aaac' + reverse_complement(ng_spacer_sequence_edit)

									else:
										spacer_oligo_top = 'caccG' + ng_spacer_sequence_edit
										spacer_oligo_bottom = 'aaac' + reverse_complement('G' + ng_spacer_sequence_edit)

									if filter_homopolymer_ts:

										if 'TTTT' not in ng_spacer_sequence_edit:
											f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + [''] * len(SILENT_CODON_COLUMNS) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')

									else:
										f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + [''] * len(SILENT_CODON_COLUMNS) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')

							counter += 1

							if number_of_pegrnas and peg_count == number_of_pegrnas:
								break

				else:

					for pegid in list(pe_design[target_name].keys())[:number_of_pegrnas or None]:

						# Write pegRNAs
						for pegRNA_entry in pe_design[target_name][pegid][0]:
							pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, pe_strand, pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries = pegRNA_entry

							pegRNA_ext_first_base = pegRNA_ext[0]
							spacer_gc_content = gc_content(pe_spacer_sequence)
							pbs_gc_content = gc_content(pegRNA_ext[rtt_length:])
							rtt_gc_content = gc_content(pegRNA_ext[:rtt_length])

							if pe_spacer_sequence[0].upper() == 'G':
								spacer_oligo_top = 'cacc' + pe_spacer_sequence + 'gtttt'
								spacer_oligo_bottom = 'ctctaaaac' + reverse_complement(pe_spacer_sequence)

							else:
								spacer_oligo_top = 'caccG' + pe_spacer_sequence + 'gtttt'
								spacer_oligo_bottom = 'ctctaaaac' + reverse_complement('G' + pe_spacer_sequence)

							pegext_oligo_top = 'gtgc' + pegRNA_ext
							pegext_oligo_bottom = 'aaaa' + reverse_complement(pegRNA_ext)

							if filter_homopolymer_ts:

								if 'TTTT' not in pe_spacer_sequence:
									f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + codon_columns(target_design[target_name], silent_mutation_selected) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')

							else:

								f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + codon_columns(target_design[target_name], silent_mutation_selected) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')

							# Sort ngRNAs

							# Write ngRNAs
							for ngRNA_entry in ngRNA_entries[:number_of_ngrnas]:
								ng_code, ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, ng_strand, nick_distance = ngRNA_entry

								spacer_gc_content = gc_content(ng_spacer_sequence_edit)

								if ng_spacer_sequence_edit[0].upper() == 'G':
									spacer_oligo_top = 'cacc' + ng_spacer_sequence_edit
									spacer_oligo_bottom = 'aaac' + reverse_complement(ng_spacer_sequence_edit)

								else:
									spacer_oligo_top = 'caccG' + ng_spacer_sequence_edit
									spacer_oligo_bottom = 'aaac' + reverse_complement('G' + ng_spacer_sequence_edit)

								if filter_homopolymer_ts:

									if 'TTTT' not in ng_spacer_sequence_edit:
										f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + [''] * len(SILENT_CODON_COLUMNS) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')

								else:
									f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + [''] * len(SILENT_CODON_COLUMNS) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')

							counter += 1

			else:
		
				for pegid in pe_design[target_name]:

					ng_continue = True

					# Write pegRNAs
					for pegRNA_entry in pe_design[target_name][pegid][0]:
						pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, pe_strand, pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit, silent_mutation_selected, ngRNA_entries = pegRNA_entry

						pegRNA_ext_first_base = pegRNA_ext[0]
						spacer_gc_content = gc_content(pe_spacer_sequence)
						pbs_gc_content = gc_content(pegRNA_ext[rtt_length:])
						rtt_gc_content = gc_content(pegRNA_ext[:rtt_length])

						if pe_spacer_sequence[0].upper() == 'G':
							spacer_oligo_top = 'cacc' + pe_spacer_sequence + 'gtttt'
							spacer_oligo_bottom = 'ctctaaaac' + reverse_complement(pe_spacer_sequence)

						else:
							spacer_oligo_top = 'caccG' + pe_spacer_sequence + 'gtttt'
							spacer_oligo_bottom = 'ctctaaaac' + reverse_complement('G' + pe_spacer_sequence)

						pegext_oligo_top = 'gtgc' + pegRNA_ext
						pegext_oligo_bottom = 'aaaa' + reverse_complement(pegRNA_ext)

						if filter_homopolymer_ts:

							if 'TTTT' not in pe_spacer_sequence:
								f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + codon_columns(target_design[target_name], silent_mutation_selected) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')

							else:
								ng_continue = False

						else:

							f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + codon_columns(target_design[target_name], silent_mutation_selected) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')

					# Write ngRNAs
						if ng_continue:
							for ngRNA_entry in ngRNA_entries[:number_of_ngrnas]:
								ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, ng_strand, nick_distance = ngRNA_entry

								spacer_gc_content = gc_content(ng_spacer_sequence_edit)

								if ng_spacer_sequence_edit[0].upper() == 'G':
									spacer_oligo_top = 'cacc' + ng_spacer_sequence_edit
									spacer_oligo_bottom = 'aaac' + reverse_complement(ng_spacer_sequence_edit)

								else:
									spacer_oligo_top = 'caccG' + ng_spacer_sequence_edit
									spacer_oligo_bottom = 'aaac' + reverse_complement('G' + ng_spacer_sequence_edit)

								if filter_homopolymer_ts:

									if 'TTTT' not in ng_spacer_sequence_edit:
										f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + [''] * len(SILENT_CODON_COLUMNS) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')

								else:

									f.write(','.join(map(str, [target_design[target_name].get('output_target_name', target_name), target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit] + [''] * len(SILENT_CODON_COLUMNS) + nucleotide_change_columns(reference_sequence, edit_sequence))) + '\n')

						# one number per pegRNA row, not per spacer: a spacer carries several silent
						# mutation alleles, and an ngRNA belongs to the allele it was cut from
						counter += 1


if __name__ == '__main__':
	main()
