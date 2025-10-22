# Design pegRNAs and ngRNAs for prime editing
##### Import libraries
import os
import sys
import re
import time
import argparse
import logging
from argparse import RawTextHelpFormatter

##### Argument handeling
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
parser.add_argument('-sat_mut', '--saturation_mutagenesis', default = False, choices = ['aa', 'aa_subs', 'aa_ins', 'aa_dels', 'base'], type = str, help = '***** Saturation mutagenesis design with prime editing. The \'aa\' option makes all amino acid substitutions (\'aa_subs\'),  +1 amino acid insertions (\'aa_ins\'), and -1 amino acid deletions (\'aa_dels\'). The \'base\' option makes DNA base changes. (Default: False) *****\n\n')
parser.add_argument('-n_pegrnas', '--number_of_pegrnas', default = 3, type = int, help = '***** The maximum number of pegRNAs to design for each input sequence. The pegRNAs are ranked by 1) PAM disrupted > PAM intact then 2) distance to edit. (Default: 3) *****\n\n')
parser.add_argument('-n_ngrnas', '--number_of_ngrnas', default = 3, type = int, help = '***** The maximum number of ngRNAs to design for each input sequence. The ngRNAs are ranked by 1) PE3b-seed > PE3b-nonseed > PE3 then 2) deviation from nicking_distance_pooled. (Default: 3) *****\n\n')
parser.add_argument('-nick_dist_pooled', '--nicking_distance_pooled', default = 75, type = int, help = '***** The nicking distance between pegRNAs and ngRNAs for pooled designs. PE3b annotation is priority, followed by nicking distance closest to this parameter. (Default: 75 bp) *****\n\n')
parser.add_argument('-homology_downstream', '--homology_downstream', default = 15, type = int, help = '***** This parameter determines the minimum RT extension length downstream of an edit for pegRNA designs. (Default: 15) *****\n\n')
parser.add_argument('-pbs_pooled', '--pbs_length_pooled', type = int, default = 14, help = '***** The PBS length to design pegRNAs for pooled design applications. (Default: 14 nt) *****\n\n')
parser.add_argument('-rtt_pooled', '--rtt_max_length_pooled', type = int, default = 50, help = '***** The maximum RTT length to design pegRNAs for pooled design applications. (Default: 50 nt) *****\n\n')


# Output directory
parser.add_argument('-out', '--out_dir', default = '0', type = str, help = '***** Name of output directory (Default: ./DATETIMESTAMP_PrimeDesign). *****\n\n')

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
if out_dir == '0':
	out_dir = '%s_PrimeDesign' % str(time.strftime("%y%m%d_%H.%M.%S", time.localtime()))

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
	for edit_idx in edit_idxs:
		edit = input_sequence[edit_idx[0]:edit_idx[1]]

		# Create edit format and number to sequence map
		if '/' in edit:
			editformat2sequence[edit] = [edit.split('/')[0].replace('(',''), edit.split('/')[1].replace(')','').lower(), edit_counter]
			editnumber2sequence[edit_counter] = [edit.split('/')[0].replace('(',''), edit.split('/')[1].replace(')','').lower()]

		elif '+' in edit:
			editformat2sequence[edit] = ['' , edit.split('+')[1].replace(')','').lower(), edit_counter]
			editnumber2sequence[edit_counter] = ['' , edit.split('+')[1].replace(')','').lower()]

		elif '-' in edit:
			editformat2sequence[edit] = [edit.split('-')[1].replace(')',''), '', edit_counter]
			editnumber2sequence[edit_counter] = [edit.split('-')[1].replace(')',''), '']

		edit_counter += 1

	edit_start = min([i.start() for i in re.finditer(r'\(', input_sequence)])
	edit_stop = max([i.start() for i in re.finditer(r'\)', input_sequence)])

	edit_span_sequence_w_ref = input_sequence[edit_start:edit_stop + 1]
	edit_span_sequence_w_edit = input_sequence[edit_start:edit_stop + 1]
	for edit in editformat2sequence:
		edit_span_sequence_w_ref = edit_span_sequence_w_ref.replace(edit, editformat2sequence[edit][0])
		edit_span_sequence_w_edit = edit_span_sequence_w_edit.replace(edit, editformat2sequence[edit][1])

	edit_start_in_ref = re.search(r'\(', input_sequence).start()
	edit_stop_in_ref_rev = re.search(r'\)', input_sequence[::-1]).start()

	edit_span_length_w_ref = len(edit_span_sequence_w_ref)
	edit_span_length_w_edit = len(edit_span_sequence_w_edit)

	reference_sequence = input_sequence
	edit_sequence = input_sequence
	editnumber_sequence = input_sequence
	for edit in editformat2sequence:
		reference_sequence = reference_sequence.replace(edit, editformat2sequence[edit][0])
		edit_sequence = edit_sequence.replace(edit, editformat2sequence[edit][1])
		editnumber_sequence = editnumber_sequence.replace(edit, str(editformat2sequence[edit][2]))

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

	if (sm_type == 'aa') | (sm_type == 'aa_subs'): # All Substitutions

		for base_index in range(0, len(sequence_to_edit), 3):

			codon_ref = sequence_to_edit[base_index:base_index + 3]

			if len(codon_ref) == 3:

				aa_ref = codon_dict[codon_ref][1]
				inner_sequence_left = sequence_to_edit[:base_index]
				inner_sequence_right = sequence_to_edit[base_index + 3:]

				aa_edit_list = [x for x in aa2codon if x != aa_ref]
				for aa_edit in aa_edit_list:

					codon_edit = aa2codon[aa_edit][0][0]
					sm_target_name_list.append('%s_%s_%sto%s' % (target_name, str(int(base_index/3 + 1)), aa_ref, aa_edit))
					sm_target_sequence_list.append(sequence_left + inner_sequence_left + '(%s/%s)' % (codon_ref, codon_edit) + inner_sequence_right + sequence_right)
	
	if (sm_type == 'aa') | (sm_type == 'aa_ins'): # All +1 Insertions
		
		for base_index in range(0, len(sequence_to_edit), 3):

			codon_ref = sequence_to_edit[base_index:base_index + 3]

			if len(codon_ref) == 3:

				aa_ref = codon_dict[codon_ref][1]
				inner_sequence_left = sequence_to_edit[:base_index]
				inner_sequence_right = sequence_to_edit[base_index+3:]

				aa_edit_list = [x for x in aa2codon if x != '*']
				for aa_edit in aa_edit_list:

					codon_edit = aa2codon[aa_edit][0][0]
					sm_target_name_list.append('%s_%s_%sto%s%s' % (target_name, str(int(base_index/3 + 1)), aa_ref, aa_ref, aa_edit))
					sm_target_sequence_list.append(sequence_left + inner_sequence_left + codon_ref + '(+%s)' % (codon_edit) + inner_sequence_right + sequence_right)
	
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

				sm_target_name_list.append('%s_%s_%s%sto%s' % (target_name, str(int(base_index/3 + 1)), aa_ref, aa2_ref, aa2_ref))
				sm_target_sequence_list.append(sequence_left + inner_sequence_left + '(-%s)' % (codon_ref) + inner_sequence_right + sequence_right)

	if sm_type == 'base': # All Base Substitutions

		base_list = ['A','T','C','G']
		for base_index in range(len(sequence_to_edit)):

			base_ref = sequence_to_edit[base_index]
			inner_sequence_left = sequence_to_edit[:base_index]
			inner_sequence_right = sequence_to_edit[base_index + 1:]

			base_edit_list = [x for x in base_list if x != base_ref.upper()]
			for base_edit in base_edit_list:

				sm_target_name_list.append('%s_%s_%sto%s' % (target_name, str(base_index+1), base_ref, base_edit))
				sm_target_sequence_list.append(sequence_left + inner_sequence_left + '(%s/%s)' % (base_ref, base_edit) + inner_sequence_right + sequence_right)

	return(sm_target_name_list, sm_target_sequence_list)

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

		target_sequence = target_sequence.upper()
		
		if saturation_mutagenesis:

			sm_target_name_list, sm_target_sequence_list = saturating_mutagenesis_input_sequences(target_name, target_sequence, saturation_mutagenesis)
			
			for sm_target_name, sm_target_sequence in zip(sm_target_name_list, sm_target_sequence_list):

				editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev = process_sequence(sm_target_sequence)

				# Initialize dictionary for the design of pegRNA spacers for each target sequence and intended edit(s)
				target_design[sm_target_name] = {'target_sequence':sm_target_sequence, 'editformat2sequence': editformat2sequence, 'editnumber2sequence': editnumber2sequence, 'reference_sequence': reference_sequence, 'edit_sequence': edit_sequence, 'editnumber_sequence': editnumber_sequence, 'edit_span_length': [edit_span_length_w_ref, edit_span_length_w_edit], 'edit_start_in_ref': edit_start_in_ref, 'edit_stop_in_ref_rev': edit_stop_in_ref_rev, 'pegRNA':{'+':[], '-':[]}, 'ngRNA':{'+':[], '-':[]}}

		else:
			editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev = process_sequence(target_sequence)

			# Initialize dictionary for the design of pegRNA spacers for each target sequence and intended edit(s)
			target_design[target_name] = {'target_sequence':target_sequence, 'editformat2sequence': editformat2sequence, 'editnumber2sequence': editnumber2sequence, 'reference_sequence': reference_sequence, 'edit_sequence': edit_sequence, 'editnumber_sequence': editnumber_sequence, 'edit_span_length': [edit_span_length_w_ref, edit_span_length_w_edit], 'edit_start_in_ref': edit_start_in_ref, 'edit_stop_in_ref_rev': edit_stop_in_ref_rev, 'pegRNA':{'+':[], '-':[]}, 'ngRNA':{'+':[], '-':[]}}

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
					silent_mutation_relative_to_edit = None
					if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]'): # Generate silent mutations

						if pe_annotate_constant == 'PAM_intact': # Generate silent mutation to distrupt PAM sequence
							
							if nick_aa_index == 0:
								codon_start_idx = 3
								codon_end_idx = 6
								original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
								aa_identity = codon_dict[original_codon][1]

								for codon_substitute in aa2codon[aa_identity]:
									if not re.search(pam_search, codon_substitute[0]):

										new_codon = codon_substitute[0].lower()
										pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
										pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_max_length_pooled])
										pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + new_codon
										pe_annotate = 'PAM_disrupted_with_silent_mutation'
										silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
										if nick2lastedit_length <= codon_start_idx: silent_mutation_relative_to_edit = 'downstream'
										elif nick2edit_length >= codon_end_idx: silent_mutation_relative_to_edit = 'upstream'
										else: silent_mutation_relative_to_edit = 'overlap'
										break

								if pe_annotate != 'PAM_disrupted_with_silent_mutation': # Generate silent mutation that does not disrupt PAM sequence
									
									for codon_substitute in aa2codon[aa_identity]:
										if codon_substitute[0] != original_codon:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
											pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_max_length_pooled])
											pe_silent_mutation = original_codon + '-to-' + new_codon
											pe_annotate = 'PAM_intact_with_silent_mutation'
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
											if nick2lastedit_length <= codon_start_idx: silent_mutation_relative_to_edit = 'downstream'
											elif nick2edit_length >= codon_end_idx: silent_mutation_relative_to_edit = 'upstream'
											else: silent_mutation_relative_to_edit = 'overlap'
											break

							elif nick_aa_index == 1:
								codon_start_idx_1 = 2
								codon_end_idx_1 = 5
								codon_start_idx_2 = 5
								codon_end_idx_2 = 8
								
								original_codon_1 = edit_sequence[pe_nick_edit_idx + codon_start_idx_1:pe_nick_edit_idx + codon_end_idx_1].upper()
								original_codon_2 = edit_sequence[pe_nick_edit_idx + codon_start_idx_2:pe_nick_edit_idx + codon_end_idx_2].upper()

								aa_identity_1 = codon_dict[original_codon_1][1]
								aa_identity_2 = codon_dict[original_codon_2][1]
								
								for codon_substitute in aa2codon[aa_identity_1]: # Try to introduce silent mutation in codon 1 first

									new_codons = codon_substitute[0] + original_codon_2
									pam_slice = new_codons[1:4]

									if not re.search(pam_search, pam_slice):

										new_codon = codon_substitute[0].lower()
										pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + rtt_length])
										pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + rtt_max_length_pooled])
										pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + pam_slice.lower()
										pe_annotate = 'PAM_disrupted_with_silent_mutation'
										silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:]
										if nick2lastedit_length <= codon_start_idx_1: silent_mutation_relative_to_edit = 'downstream'
										elif nick2edit_length >= codon_end_idx_1: silent_mutation_relative_to_edit = 'upstream'
										else: silent_mutation_relative_to_edit = 'overlap'
										break

								if pe_annotate != 'PAM_disrupted_with_silent_mutation': # Try to introduce silent mutation in codon 2 if codon 1 did not work
									
									for codon_substitute in aa2codon[aa_identity_2]:

										new_codons = original_codon_1 + codon_substitute[0]
										pam_slice = new_codons[1:4]

										if not re.search(pam_search, pam_slice):

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + rtt_length])
											pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + rtt_max_length_pooled])
											pe_silent_mutation = pe_pam_ref + '-to-' + pam_slice.lower()
											pe_annotate = 'PAM_disrupted_with_silent_mutation'
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:]
											if nick2lastedit_length <= codon_start_idx_2: silent_mutation_relative_to_edit = 'downstream'
											elif nick2edit_length >= codon_end_idx_2: silent_mutation_relative_to_edit = 'upstream'
											else: silent_mutation_relative_to_edit = 'overlap'
											break
								
								if 'silent_mutation' not in pe_annotate: # Generate silent mutation that does not disrupt PAM sequence
									
									for codon_substitute in aa2codon[aa_identity_1]: # Try to introduce silent mutation in codon 1 first
										if codon_substitute[0] != original_codon_1:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + rtt_length])
											pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + rtt_max_length_pooled])
											pe_silent_mutation = original_codon_1 + '-to-' + new_codon
											pe_annotate = 'PAM_intact_with_silent_mutation'
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:]
											if nick2lastedit_length <= codon_start_idx_1: silent_mutation_relative_to_edit = 'downstream'
											elif nick2edit_length >= codon_end_idx_1: silent_mutation_relative_to_edit = 'upstream'
											else: silent_mutation_relative_to_edit = 'overlap'
											break
								
								if 'silent_mutation' not in pe_annotate: # Try to introduce silent mutation in codon 2 if codon 1 did not work

									for codon_substitute in aa2codon[aa_identity_2]:
										if codon_substitute[0] != original_codon_2:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + rtt_length])
											pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + rtt_max_length_pooled])
											pe_silent_mutation = original_codon_2 + '-to-' + new_codon
											pe_annotate = 'PAM_intact_with_silent_mutation'
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:]
											if nick2lastedit_length <= codon_start_idx_2: silent_mutation_relative_to_edit = 'downstream'
											elif nick2edit_length >= codon_end_idx_2: silent_mutation_relative_to_edit = 'upstream'
											else: silent_mutation_relative_to_edit = 'overlap'
											break

							elif nick_aa_index == 2:
								codon_start_idx = 4
								codon_end_idx = 7
								original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
								aa_identity = codon_dict[original_codon][1]

								for codon_substitute in aa2codon[aa_identity]:

									pam_slice = edit_sequence[pe_nick_edit_idx + 3:pe_nick_edit_idx + 4].upper() + codon_substitute[0][:2]
									if not re.search(pam_search, pam_slice):

										new_codon = codon_substitute[0].lower()
										pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
										pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_max_length_pooled])
										pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + pam_slice.lower()
										pe_annotate = 'PAM_disrupted_with_silent_mutation'
										silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
										if nick2lastedit_length <= codon_start_idx: silent_mutation_relative_to_edit = 'downstream'
										elif nick2edit_length >= codon_end_idx: silent_mutation_relative_to_edit = 'upstream'
										else: silent_mutation_relative_to_edit = 'overlap'
										break

								if pe_annotate != 'PAM_disrupted_with_silent_mutation': # Generate silent mutation that does not disrupt PAM sequence
									
									for codon_substitute in aa2codon[aa_identity]:
										if codon_substitute[0] != original_codon:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
											pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_max_length_pooled])
											pe_silent_mutation = original_codon + '-to-' + new_codon
											pe_annotate = 'PAM_intact_with_silent_mutation'
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
											if nick2lastedit_length <= codon_start_idx: silent_mutation_relative_to_edit = 'downstream'
											elif nick2edit_length >= codon_end_idx: silent_mutation_relative_to_edit = 'upstream'
											else: silent_mutation_relative_to_edit = 'overlap'
											break
						
						if pe_annotate_constant == 'PAM_disrupted' or silent_mutation_relative_to_edit == 'overlap': # Generate silent mutation 5' or 3' of the PAM

							if nick_aa_index == 0:
								codon_start_idx = 0 # Try to introduce silent mutation 5' of the PAM first
								codon_end_idx = 3
								original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
								aa_identity = codon_dict[original_codon][1]

								if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
									for codon_substitute in aa2codon[aa_identity]:
										if codon_substitute[0] != original_codon:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
											pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_max_length_pooled])
											pe_silent_mutation = original_codon + '-to-' + new_codon

											if silent_mutation_relative_to_edit == 'overlap': 
												if 'PAM_disrupted' in pe_annotate: 
													pe_annotate = 'silent_mutation_and_PAM_disrupted'
												else:
													pe_annotate = 'silent_mutation_and_PAM_intact'
											else:
												pe_annotate = 'silent_mutation_and_PAM_disrupted'
											
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
											silent_mutation_relative_to_edit = 'upstream'
											break

								else: # If methionine or tryptophan codons, try to introduce silent mutation 3' of the PAM
									codon_start_idx = 6
									codon_end_idx = 9
									original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
									aa_identity = codon_dict[original_codon][1]

									if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
										for codon_substitute in aa2codon[aa_identity]:
											if codon_substitute[0] != original_codon:
												
												new_codon = codon_substitute[0].lower()
												pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
												pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_max_length_pooled])
												pe_silent_mutation = original_codon + '-to-' + new_codon

												if silent_mutation_relative_to_edit == 'overlap': 
													if 'PAM_disrupted' in pe_annotate: 
														pe_annotate = 'PAM_disrupted_and_silent_mutation'
													else:
														pe_annotate = 'PAM_intact_and_silent_mutation'
												else:
													pe_annotate = 'PAM_disrupted_and_silent_mutation'
												
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
												silent_mutation_relative_to_edit = 'downstream'
												break

							elif nick_aa_index == 1:
								codon_start_idx = -1 # Try to introduce silent mutation 5' of the PAM first
								codon_end_idx = 2
								original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
								aa_identity = codon_dict[original_codon][1]

								if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
									for codon_substitute in aa2codon[aa_identity]:
										if (codon_substitute[0][0] == original_codon[0][0]) & (codon_substitute[0] != original_codon): # Keep the first base the same to preserve PBS sequence

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
											pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_max_length_pooled])
											pe_silent_mutation = original_codon + '-to-' + new_codon
											
											if silent_mutation_relative_to_edit == 'overlap': 
												if 'PAM_disrupted' in pe_annotate: 
													pe_annotate = 'silent_mutation_and_PAM_disrupted'
												else:
													pe_annotate = 'silent_mutation_and_PAM_intact'
											else:
												pe_annotate = 'silent_mutation_and_PAM_disrupted'

											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
											silent_mutation_relative_to_edit = 'upstream'
											break
								
								else: # If methionine or tryptophan codons, try to introduce silent mutation 3' of the PAM
									codon_start_idx = 8
									codon_end_idx = 11
									original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
									aa_identity = codon_dict[original_codon][1]

									if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
										for codon_substitute in aa2codon[aa_identity]:
											if codon_substitute[0] != original_codon:
												
												new_codon = codon_substitute[0].lower()
												pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
												pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_max_length_pooled])
												pe_silent_mutation = original_codon + '-to-' + new_codon
												
												if silent_mutation_relative_to_edit == 'overlap': 
													if 'PAM_disrupted' in pe_annotate: 
														pe_annotate = 'PAM_disrupted_and_silent_mutation'
													else:
														pe_annotate = 'PAM_intact_and_silent_mutation'
												else:
													pe_annotate = 'PAM_disrupted_and_silent_mutation'
												
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
												silent_mutation_relative_to_edit = 'downstream'
												break

							elif nick_aa_index == 2:
								codon_start_idx = 1 # Try to introduce silent mutation 5' of the PAM first
								codon_end_idx = 4
								original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
								aa_identity = codon_dict[original_codon][1]

								if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
									for codon_substitute in aa2codon[aa_identity]:
										if codon_substitute[0] != original_codon:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
											pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_max_length_pooled])
											pe_silent_mutation = original_codon + '-to-' + new_codon
											
											if silent_mutation_relative_to_edit == 'overlap': 
												if 'PAM_disrupted' in pe_annotate: 
													pe_annotate = 'silent_mutation_and_PAM_disrupted'
												else:
													pe_annotate = 'silent_mutation_and_PAM_intact'
											else:
												pe_annotate = 'silent_mutation_and_PAM_disrupted'
												
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
											silent_mutation_relative_to_edit = 'upstream'
											break
								
								else: # If methionine or tryptophan codons, try to introduce silent mutation 3' of the PAM
									codon_start_idx = 7
									codon_end_idx = 10
									original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
									aa_identity = codon_dict[original_codon][1]

									if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
										for codon_substitute in aa2codon[aa_identity]:
											if codon_substitute[0] != original_codon:

												new_codon = codon_substitute[0].lower()
												pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
												pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_max_length_pooled])
												pe_silent_mutation = original_codon + '-to-' + new_codon
												
												if silent_mutation_relative_to_edit == 'overlap': 
													if 'PAM_disrupted' in pe_annotate: 
														pe_annotate = 'PAM_disrupted_and_silent_mutation'
													else:
														pe_annotate = 'PAM_intact_and_silent_mutation'
												else:
													pe_annotate = 'PAM_disrupted_and_silent_mutation'
												
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
												silent_mutation_relative_to_edit = 'downstream'
												break
					
					# Create pegRNA ID based on annotation
					if 'PAM_disrupted' in pe_annotate and 'silent_mutation' in pe_annotate:
						pe_annotate_code = 0
					elif 'silent_mutation' in pe_annotate:
						pe_annotate_code = 1
					elif 'PAM_disrupted' in pe_annotate:
						pe_annotate_code = 2
					else:
						pe_annotate_code = 3

					pegid = '_'.join(map(str, [str(pe_annotate_code) + '0'*(3 - len(str(abs(nick2lastedit_length)))) + str(abs(nick2lastedit_length)), pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+']))
					
					# Check to see if pegRNA extension is within input sequence
					if len(pegRNA_ext) == (pbs_length + rtt_length):

						# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
						if pegid not in pe_design[target_name]:

							# First list is for peg extension, second list is for nicking guide
							pe_design[target_name][pegid] = [[],[]]
						
						# Store pegRNA design
						if silent_mutation_edit_sequence == '':
							if pe_pam_ref_silent_mutation == '':
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])

							else:
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref_silent_mutation, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])
						else:
							if pe_pam_ref_silent_mutation == '':
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit])

							else:
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref_silent_mutation, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit])

					# Create ngRNAs targeting (-) strand for (+) pegRNAs
					if pegid in pe_design[target_name]:
						for ng_minus in target_design[target_name]['ngRNA']['-']:
							ng_nick_ref_idx, ng_edit_start_idx, ng_edit_end_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_minus
							nick_distance = ng_nick_ref_idx - pe_nick_ref_idx

							if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]') and (len(silent_mutation_edit_sequence) > 0):
								ng_spacer_sequence_edit = silent_mutation_edit_sequence[ng_edit_start_idx:ng_edit_end_idx]

								mutation_indices = [i for i, a in enumerate(ng_spacer_sequence_edit) if a.islower()]
								if len(mutation_indices) > 0:
									if len([1 for x in mutation_indices if x < 10]) > 0:
										ng_annotate = 'PE3b-seed'

									else:
										ng_annotate = 'PE3b-nonseed'
								else:
									ng_annotate = 'PE3'

							if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):

								if ng_annotate == 'PE3b-seed':
									ng_code = 0
								elif ng_annotate == 'PE3b-nonseed':
									ng_code = 1
								else:
									ng_code = 2

								pe_design[target_name][pegid][1].append([str(ng_code) + '0'*(3 - len(str(abs(abs(nick_distance) - nicking_distance_pooled)))) + str(abs(abs(nick_distance) - nicking_distance_pooled)), ng_nick_ref_idx, reverse_complement(ng_spacer_sequence_edit), reverse_complement(ng_pam_edit), ng_annotate, '-', nick_distance])

						pe_design[target_name][pegid][1] = sorted(pe_design[target_name][pegid][1])

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
					silent_mutation_relative_to_edit = None
					if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]'):
						
						if pe_annotate_constant == 'PAM_intact':
							
							if nick_aa_index == 0:
								codon_start_idx = -6
								codon_end_idx = -3
								original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
								aa_identity = codon_dict[original_codon][1]

								for codon_substitute in aa2codon[aa_identity]:
									if not re.search(reverse_complement(pam_search), codon_substitute[0]):

										new_codon = codon_substitute[0].lower()
										pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
										pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
										pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(new_codon)
										pe_annotate = 'PAM_disrupted_with_silent_mutation'
										silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
										if nick2lastedit_length >= abs(codon_start_idx): silent_mutation_relative_to_edit = 'downstream'
										elif nick2edit_length <= abs(codon_end_idx): silent_mutation_relative_to_edit = 'upstream'
										else: silent_mutation_relative_to_edit = 'overlap'
										break
								
								if pe_annotate != 'PAM_disrupted_with_silent_mutation': # Generate silent mutation that does not disrupt PAM sequence
									
									for codon_substitute in aa2codon[aa_identity]:
										if codon_substitute[0] != original_codon:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
											pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
											pe_silent_mutation = original_codon + '-to-' + new_codon
											pe_annotate = 'PAM_intact_with_silent_mutation'
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
											if nick2lastedit_length >= abs(codon_start_idx): silent_mutation_relative_to_edit = 'downstream'
											elif nick2edit_length <= abs(codon_end_idx): silent_mutation_relative_to_edit = 'upstream'
											else: silent_mutation_relative_to_edit = 'overlap'
											break

							elif nick_aa_index == 1:
								codon_start_idx = -7
								codon_end_idx = -4
								original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
								aa_identity = codon_dict[original_codon][1]

								for codon_substitute in aa2codon[aa_identity]:

									pam_slice = codon_substitute[0][1:] + edit_sequence[pe_nick_edit_idx - 4:pe_nick_edit_idx - 3].upper()
									if not re.search(reverse_complement(pam_search), pam_slice):

										new_codon = codon_substitute[0].lower()
										pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
										pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
										pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(pam_slice).lower()
										pe_annotate = 'PAM_disrupted_with_silent_mutation'
										silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
										if nick2lastedit_length >= abs(codon_start_idx): silent_mutation_relative_to_edit = 'downstream'
										elif nick2edit_length <= abs(codon_end_idx): silent_mutation_relative_to_edit = 'upstream'
										else: silent_mutation_relative_to_edit = 'overlap'
										break
								
								if pe_annotate != 'PAM_disrupted_with_silent_mutation': # Generate silent mutation that does not disrupt PAM sequence
									
									for codon_substitute in aa2codon[aa_identity]:
										if codon_substitute[0] != original_codon:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
											pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
											pe_silent_mutation = original_codon + '-to-' + new_codon
											pe_annotate = 'PAM_intact_with_silent_mutation'
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
											if nick2lastedit_length >= abs(codon_start_idx): silent_mutation_relative_to_edit = 'downstream'
											elif nick2edit_length <= abs(codon_end_idx): silent_mutation_relative_to_edit = 'upstream'
											else: silent_mutation_relative_to_edit = 'overlap'
											break

							elif nick_aa_index == 2:
								codon_start_idx_1 = -8
								codon_end_idx_1 = -5
								codon_start_idx_2 = -5
								codon_end_idx_2 = -2

								original_codon_1 = edit_sequence[pe_nick_edit_idx + codon_start_idx_1:pe_nick_edit_idx + codon_end_idx_1].upper()
								original_codon_2 = edit_sequence[pe_nick_edit_idx + codon_start_idx_2:pe_nick_edit_idx + codon_end_idx_2].upper()

								aa_identity_1 = codon_dict[original_codon_1][1]
								aa_identity_2 = codon_dict[original_codon_2][1]

								for codon_substitute in aa2codon[aa_identity_1]:

									new_codons = codon_substitute[0] + original_codon_2
									pam_slice = new_codons[2:5]

									if not re.search(reverse_complement(pam_search), pam_slice):

										new_codon = codon_substitute[0].lower()
										pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + pbs_length]
										pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + pbs_length]
										pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(pam_slice).lower()
										pe_annotate = 'PAM_disrupted_with_silent_mutation'
										silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:]
										if nick2lastedit_length >= abs(codon_start_idx_1): silent_mutation_relative_to_edit = 'downstream'
										elif nick2edit_length <= abs(codon_end_idx_1): silent_mutation_relative_to_edit = 'upstream'
										else: silent_mutation_relative_to_edit = 'overlap'
										break

								if pe_annotate != 'PAM_disrupted_with_silent_mutation':

									for codon_substitute in aa2codon[aa_identity_2]:

										new_codons = original_codon_1 + codon_substitute[0]
										pam_slice = new_codons[2:5]

										if not re.search(reverse_complement(pam_search), pam_slice):

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + pbs_length]
											pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + pbs_length]
											pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(pam_slice).lower()
											pe_annotate = 'PAM_disrupted_with_silent_mutation'
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:]
											if nick2lastedit_length >= abs(codon_start_idx_2): silent_mutation_relative_to_edit = 'downstream'
											elif nick2edit_length <= abs(codon_end_idx_2): silent_mutation_relative_to_edit = 'upstream'
											else: silent_mutation_relative_to_edit = 'overlap'
											break
								
								if 'silent_mutation' not in pe_annotate: # Generate silent mutation that does not disrupt PAM sequence
									
									for codon_substitute in aa2codon[aa_identity_1]: # Try to introduce silent mutation in codon 1 first
										if codon_substitute[0] != original_codon_1:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + pbs_length]
											pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + pbs_length]
											pe_silent_mutation = original_codon_1 + '-to-' + new_codon
											pe_annotate = 'PAM_intact_with_silent_mutation'
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:]
											if nick2lastedit_length >= abs(codon_start_idx_1): silent_mutation_relative_to_edit = 'downstream'
											elif nick2edit_length <= abs(codon_end_idx_1): silent_mutation_relative_to_edit = 'upstream'
											else: silent_mutation_relative_to_edit = 'overlap'
											break
									
								if 'silent_mutation' not in pe_annotate: # Try to introduce silent mutation in codon 2 if codon 1 did not work

									for codon_substitute in aa2codon[aa_identity_2]:
										if codon_substitute[0] != original_codon_2:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + pbs_length]
											pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + pbs_length]
											pe_silent_mutation = original_codon_2 + '-to-' + new_codon
											pe_annotate = 'PAM_intact_with_silent_mutation'
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:]
											if nick2lastedit_length >= abs(codon_start_idx_2): silent_mutation_relative_to_edit = 'downstream'
											elif nick2edit_length <= abs(codon_end_idx_2): silent_mutation_relative_to_edit = 'upstream'
											else: silent_mutation_relative_to_edit = 'overlap'
											break

						if pe_annotate_constant == 'PAM_disrupted' or silent_mutation_relative_to_edit == 'overlap': # Generate silent mutation 5' or 3' of the PAM
							
							if nick_aa_index == 0:
								codon_start_idx = -3 # Try to introduce silent mutation 3' of the PAM first
								codon_end_idx = 0
								original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
								aa_identity = codon_dict[original_codon][1]

								if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
									for codon_substitute in aa2codon[aa_identity]:
										if codon_substitute[0] != original_codon:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
											pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
											pe_silent_mutation = original_codon + '-to-' + new_codon
											
											if silent_mutation_relative_to_edit == 'overlap': 
												if 'PAM_disrupted' in pe_annotate: 
													pe_annotate = 'PAM_disrupted_and_silent_mutation'
												else:
													pe_annotate = 'PAM_intact_and_silent_mutation'
											else:
												pe_annotate = 'PAM_disrupted_and_silent_mutation'
											
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
											silent_mutation_relative_to_edit = 'downstream'
											break

								else: # If methionine or tryptophan codons, try to introduce silent mutation 5' of the PAM
									codon_start_idx = -9
									codon_end_idx = -6
									original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
									aa_identity = codon_dict[original_codon][1]

									if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
										for codon_substitute in aa2codon[aa_identity]:
											if codon_substitute[0] != original_codon:
												
												new_codon = codon_substitute[0].lower()
												pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
												pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
												pe_silent_mutation = original_codon + '-to-' + new_codon
												
												if silent_mutation_relative_to_edit == 'overlap': 
													if 'PAM_disrupted' in pe_annotate: 
														pe_annotate = 'silent_mutation_and_PAM_disrupted'
													else:
														pe_annotate = 'silent_mutation_and_PAM_intact'
												else:
													pe_annotate = 'silent_mutation_and_PAM_disrupted'
												
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
												silent_mutation_relative_to_edit = 'upstream'
												break
							
							elif nick_aa_index == 1:
								codon_start_idx = -4 # Try to introduce silent mutation 3' of the PAM first
								codon_end_idx = -1
								original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
								aa_identity = codon_dict[original_codon][1]

								if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
									for codon_substitute in aa2codon[aa_identity]:
										if codon_substitute[0] != original_codon:

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
											pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
											pe_silent_mutation = original_codon + '-to-' + new_codon
											
											if silent_mutation_relative_to_edit == 'overlap': 
												if 'PAM_disrupted' in pe_annotate: 
													pe_annotate = 'PAM_disrupted_and_silent_mutation'
												else:
													pe_annotate = 'PAM_intact_and_silent_mutation'
											else:
												pe_annotate = 'PAM_disrupted_and_silent_mutation'

											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
											silent_mutation_relative_to_edit = 'downstream'
											break

								else: # If methionine or tryptophan codons, try to introduce silent mutation 5' of the PAM
									codon_start_idx = -10
									codon_end_idx = -7
									original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
									aa_identity = codon_dict[original_codon][1]

									if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
										for codon_substitute in aa2codon[aa_identity]:
											if codon_substitute[0] != original_codon:
												
												new_codon = codon_substitute[0].lower()
												pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
												pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
												pe_silent_mutation = original_codon + '-to-' + new_codon
												
												if silent_mutation_relative_to_edit == 'overlap': 
													if 'PAM_disrupted' in pe_annotate: 
														pe_annotate = 'silent_mutation_and_PAM_disrupted'
													else:
														pe_annotate = 'silent_mutation_and_PAM_intact'
												else:
													pe_annotate = 'silent_mutation_and_PAM_disrupted'
												
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
												silent_mutation_relative_to_edit = 'upstream'
												break
							
							elif nick_aa_index == 2:
								codon_start_idx = -2 # Try to introduce silent mutation 3' of the PAM first
								codon_end_idx = 1
								original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
								aa_identity = codon_dict[original_codon][1]

								if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
									for codon_substitute in aa2codon[aa_identity]:
										if (codon_substitute[0][2] == original_codon[2]) & (codon_substitute[0] != original_codon): # Keep the third base the same to preserve PBS sequence

											new_codon = codon_substitute[0].lower()
											pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
											pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
											pe_silent_mutation = original_codon + '-to-' + new_codon
											
											if silent_mutation_relative_to_edit == 'overlap': 
												if 'PAM_disrupted' in pe_annotate: 
													pe_annotate = 'PAM_disrupted_and_silent_mutation'
												else:
													pe_annotate = 'PAM_intact_and_silent_mutation'
											else:
												pe_annotate = 'PAM_disrupted_and_silent_mutation'
											
											silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
											silent_mutation_relative_to_edit = 'downstream'
											break

								else: # If methionine or tryptophan codons, try to introduce silent mutation 5' of the PAM
									codon_start_idx = -11
									codon_end_idx = -8
									original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
									aa_identity = codon_dict[original_codon][1]

									if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
										for codon_substitute in aa2codon[aa_identity]:
											if codon_substitute[0] != original_codon:
												
												new_codon = codon_substitute[0].lower()
												pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
												pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
												pe_silent_mutation = original_codon + '-to-' + new_codon
												
												if silent_mutation_relative_to_edit == 'overlap': 
													if 'PAM_disrupted' in pe_annotate: 
														pe_annotate = 'silent_mutation_and_PAM_disrupted'
													else:
														pe_annotate = 'silent_mutation_and_PAM_intact'
												else:
													pe_annotate = 'silent_mutation_and_PAM_disrupted'
												
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
												silent_mutation_relative_to_edit = 'upstream'
												break

					# Create pegRNA ID based on annotation
					if 'PAM_disrupted' in pe_annotate and 'silent_mutation' in pe_annotate:
						pe_annotate_code = 0
					elif 'silent_mutation' in pe_annotate:
						pe_annotate_code = 1
					elif 'PAM_disrupted' in pe_annotate:
						pe_annotate_code = 2
					else:
						pe_annotate_code = 3

					pegid = '_'.join(map(str, [str(pe_annotate_code) + '0'*(3 - len(str(abs(nick2lastedit_length)))) + str(abs(nick2lastedit_length)), pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-']))
					
					# Check to see if pegRNA extension is within input sequence
					if len(pegRNA_ext) == (pbs_length + rtt_length):

						# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
						if pegid not in pe_design[target_name]:

							# First list is for peg extension, second list is for nicking guide
							pe_design[target_name][pegid] = [[],[]]
						
						# Store pegRNA design
						if silent_mutation_edit_sequence == '':
							if pe_pam_ref_silent_mutation == '':
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), reverse_complement(pe_pam_ref), pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])
							
							else:
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), pe_pam_ref_silent_mutation, pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])
						
						else:
							if pe_pam_ref_silent_mutation == '':
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), reverse_complement(pe_pam_ref), pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit])
							
							else:
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), pe_pam_ref_silent_mutation, pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit])

					# Create ngRNAs targeting (+) strand for (-) pegRNAs
					if pegid in pe_design[target_name]:
						for ng_plus in target_design[target_name]['ngRNA']['+']:
							ng_nick_ref_idx, ng_edit_start_idx, ng_edit_end_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_plus
							nick_distance = ng_nick_ref_idx - pe_nick_ref_idx

							if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]') and (len(silent_mutation_edit_sequence) > 0):
								ng_spacer_sequence_edit = silent_mutation_edit_sequence[ng_edit_start_idx:ng_edit_end_idx]

								mutation_indices = [i for i, a in enumerate(ng_spacer_sequence_edit) if a.islower()]
								if len(mutation_indices) > 0:
									if len([1 for x in mutation_indices if x >= 10]) > 0:
										ng_annotate = 'PE3b-seed'

									else:
										ng_annotate = 'PE3b-nonseed'
								else:
									ng_annotate = 'PE3'

							if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):

								if ng_annotate == 'PE3b-seed':
									ng_code = 0
								elif ng_annotate == 'PE3b-nonseed':
									ng_code = 1
								else:
									ng_code = 2

								pe_design[target_name][pegid][1].append([str(ng_code) + '0'*(3 - len(str(abs(abs(nick_distance) - nicking_distance_pooled)))) + str(abs(abs(nick_distance) - nicking_distance_pooled)), ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, '+', nick_distance])

						pe_design[target_name][pegid][1] = sorted(pe_design[target_name][pegid][1])

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
							silent_mutation_relative_to_edit = None
							if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]'):

								if pe_annotate_constant == 'PAM_intact':
									
									if nick_aa_index == 0:
										codon_start_idx = 3
										codon_end_idx = 6
										original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
										aa_identity = codon_dict[original_codon][1]

										for codon_substitute in aa2codon[aa_identity]:
											if not re.search(pam_search, codon_substitute[0]):

												new_codon = codon_substitute[0].lower()
												pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
												pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + new_codon
												pe_annotate = 'PAM_disrupted_with_silent_mutation'
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
												if nick2lastedit_length <= codon_start_idx: silent_mutation_relative_to_edit = 'downstream'
												elif nick2edit_length >= codon_end_idx: silent_mutation_relative_to_edit = 'upstream'
												else: silent_mutation_relative_to_edit = 'overlap'
												break
										
										if pe_annotate != 'PAM_disrupted_with_silent_mutation': # Generate silent mutation that does not disrupt PAM sequence
									
											for codon_substitute in aa2codon[aa_identity]:
												if codon_substitute[0] != original_codon:

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
													pe_silent_mutation = original_codon + '-to-' + new_codon
													pe_annotate = 'PAM_intact_with_silent_mutation'
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
													if nick2lastedit_length <= codon_start_idx: silent_mutation_relative_to_edit = 'downstream'
													elif nick2edit_length >= codon_end_idx: silent_mutation_relative_to_edit = 'upstream'
													else: silent_mutation_relative_to_edit = 'overlap'
													break

									elif nick_aa_index == 1:
										codon_start_idx_1 = 2
										codon_end_idx_1 = 5
										codon_start_idx_2 = 5
										codon_end_idx_2 = 8

										original_codon_1 = edit_sequence[pe_nick_edit_idx + codon_start_idx_1:pe_nick_edit_idx + codon_end_idx_1].upper()
										original_codon_2 = edit_sequence[pe_nick_edit_idx + codon_start_idx_2:pe_nick_edit_idx + codon_end_idx_2].upper()

										aa_identity_1 = codon_dict[original_codon_1][1]
										aa_identity_2 = codon_dict[original_codon_2][1]

										for codon_substitute in aa2codon[aa_identity_1]:

											new_codons = codon_substitute[0] + original_codon_2
											pam_slice = new_codons[1:4]

											if not re.search(pam_search, pam_slice):

												new_codon = codon_substitute[0].lower()
												pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + rtt_length])
												pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + pam_slice.lower()
												pe_annotate = 'PAM_disrupted_with_silent_mutation'
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:]
												if nick2lastedit_length <= codon_start_idx_1: silent_mutation_relative_to_edit = 'downstream'
												elif nick2edit_length >= codon_end_idx_1: silent_mutation_relative_to_edit = 'upstream'
												else: silent_mutation_relative_to_edit = 'overlap'
												break

										if pe_annotate != 'PAM_disrupted_with_silent_mutation': # Try to introduce silent mutation in codon 2 if codon 1 did not work

											for codon_substitute in aa2codon[aa_identity_2]:

												new_codons = original_codon_1 + codon_substitute[0]
												pam_slice = new_codons[1:4]

												if not re.search(pam_search, pam_slice):

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + rtt_length])
													pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + pam_slice.lower()
													pe_annotate = 'PAM_disrupted_with_silent_mutation'
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:]
													if nick2lastedit_length <= codon_start_idx_2: silent_mutation_relative_to_edit = 'downstream'
													elif nick2edit_length >= codon_end_idx_2: silent_mutation_relative_to_edit = 'upstream'
													else: silent_mutation_relative_to_edit = 'overlap'
													break

										if 'silent_mutation' not in pe_annotate: # Generate silent mutation that does not disrupt PAM sequence
											
											for codon_substitute in aa2codon[aa_identity_1]: # Try to introduce silent mutation in codon 1 first
												if codon_substitute[0] != original_codon_1:

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + rtt_length])
													pe_silent_mutation = original_codon_1 + '-to-' + new_codon
													pe_annotate = 'PAM_intact_with_silent_mutation'
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:]
													if nick2lastedit_length <= codon_start_idx_1: silent_mutation_relative_to_edit = 'downstream'
													elif nick2edit_length >= codon_end_idx_1: silent_mutation_relative_to_edit = 'upstream'
													else: silent_mutation_relative_to_edit = 'overlap'
													break
										
										if 'silent_mutation' not in pe_annotate: # Try to introduce silent mutation in codon 2 if codon 1 did not work

											for codon_substitute in aa2codon[aa_identity_2]:
												if codon_substitute[0] != original_codon_2:

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + rtt_length])
													pe_silent_mutation = original_codon_2 + '-to-' + new_codon
													pe_annotate = 'PAM_intact_with_silent_mutation'
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:]
													if nick2lastedit_length <= codon_start_idx_2: silent_mutation_relative_to_edit = 'downstream'
													elif nick2edit_length >= codon_end_idx_2: silent_mutation_relative_to_edit = 'upstream'
													else: silent_mutation_relative_to_edit = 'overlap'

													break

									elif nick_aa_index == 2:
										codon_start_idx = 4
										codon_end_idx = 7
										original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
										aa_identity = codon_dict[original_codon][1]

										for codon_substitute in aa2codon[aa_identity]:

											pam_slice = edit_sequence[pe_nick_edit_idx + 3:pe_nick_edit_idx + 4].upper() + codon_substitute[0][:2]
											if not re.search(pam_search, pam_slice):

												new_codon = codon_substitute[0].lower()
												pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
												pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + pam_slice.lower()
												pe_annotate = 'PAM_disrupted_with_silent_mutation'
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
												if nick2lastedit_length <= codon_start_idx: silent_mutation_relative_to_edit = 'downstream'
												elif nick2edit_length >= codon_end_idx: silent_mutation_relative_to_edit = 'upstream'
												else: silent_mutation_relative_to_edit = 'overlap'
												break
									
										if pe_annotate != 'PAM_disrupted_with_silent_mutation': # Generate silent mutation that does not disrupt PAM sequence
										
											for codon_substitute in aa2codon[aa_identity]:
												if codon_substitute[0] != original_codon:

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
													pe_silent_mutation = original_codon + '-to-' + new_codon
													pe_annotate = 'PAM_intact_with_silent_mutation'
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
													if nick2lastedit_length <= codon_start_idx: silent_mutation_relative_to_edit = 'downstream'
													elif nick2edit_length >= codon_end_idx: silent_mutation_relative_to_edit = 'upstream'
													else: silent_mutation_relative_to_edit = 'overlap'
													break

								if pe_annotate_constant == 'PAM_disrupted' or silent_mutation_relative_to_edit == 'overlap': # Generate silent mutation 5' or 3' of the PAM
									
									if nick_aa_index == 0:
										codon_start_idx = 0 # Try to introduce silent mutation 5' of the PAM first
										codon_end_idx = 3
										original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
										aa_identity = codon_dict[original_codon][1]

										if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons

											for codon_substitute in aa2codon[aa_identity]:
												if codon_substitute[0] != original_codon:

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
													pe_silent_mutation = original_codon + '-to-' + new_codon

													if silent_mutation_relative_to_edit == 'overlap': 
														if 'PAM_disrupted' in pe_annotate: 
															pe_annotate = 'silent_mutation_and_PAM_disrupted'
														else:
															pe_annotate = 'silent_mutation_and_PAM_intact'
													else:
														pe_annotate = 'silent_mutation_and_PAM_disrupted'
													
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
													silent_mutation_relative_to_edit = 'upstream'
													break

										else: # If methionine or tryptophan codons, try to introduce silent mutation 3' of the PAM
											codon_start_idx = 6
											codon_end_idx = 9

											original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
											aa_identity = codon_dict[original_codon][1]

											if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
												for codon_substitute in aa2codon[aa_identity]:
													if codon_substitute[0] != original_codon:
														
														new_codon = codon_substitute[0].lower()
														pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
														pe_silent_mutation = original_codon + '-to-' + new_codon

														if silent_mutation_relative_to_edit == 'overlap': 
															if 'PAM_disrupted' in pe_annotate: 
																pe_annotate = 'PAM_disrupted_and_silent_mutation'
															else:
																pe_annotate = 'PAM_intact_and_silent_mutation'
														else:
															pe_annotate = 'PAM_disrupted_and_silent_mutation'
														
														silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
														silent_mutation_relative_to_edit = 'downstream'
														break
														
											elif nick_aa_index == 1:
												codon_start_idx = -1 # Try to introduce silent mutation 5' of the PAM first
												codon_end_idx = 2
												original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
												aa_identity = codon_dict[original_codon][1]

												if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
													for codon_substitute in aa2codon[aa_identity]:
														if codon_substitute[0] != original_codon:

															new_codon = codon_substitute[0].lower()
															pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
															pe_silent_mutation = original_codon + '-to-' + new_codon
															
															if silent_mutation_relative_to_edit == 'overlap': 
																if 'PAM_disrupted' in pe_annotate: 
																	pe_annotate = 'silent_mutation_and_PAM_disrupted'
																else:
																	pe_annotate = 'silent_mutation_and_PAM_intact'
															else:
																pe_annotate = 'silent_mutation_and_PAM_disrupted'

															silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
															silent_mutation_relative_to_edit = 'upstream'
															break
												
												else: # If methionine or tryptophan codons, try to introduce silent mutation 3' of the PAM
													codon_start_idx = 8
													codon_end_idx = 11

													original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
													aa_identity = codon_dict[original_codon][1]

													if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
														for codon_substitute in aa2codon[aa_identity]:
															if codon_substitute[0] != original_codon:
																
																new_codon = codon_substitute[0].lower()
																pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
																pe_silent_mutation = original_codon + '-to-' + new_codon
																
																if silent_mutation_relative_to_edit == 'overlap': 
																	if 'PAM_disrupted' in pe_annotate: 
																		pe_annotate = 'PAM_disrupted_and_silent_mutation'
																	else:
																		pe_annotate = 'PAM_intact_and_silent_mutation'
																else:
																	pe_annotate = 'PAM_disrupted_and_silent_mutation'
																
																silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
																silent_mutation_relative_to_edit = 'downstream'
																break
											
											elif nick_aa_index == 2:
												codon_start_idx = 1 # Try to introduce silent mutation 5' of the PAM first
												codon_end_idx = 4
												original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
												aa_identity = codon_dict[original_codon][1]

												if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
													for codon_substitute in aa2codon[aa_identity]:
														if codon_substitute[0] != original_codon:

															new_codon = codon_substitute[0].lower()
															pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
															pe_silent_mutation = original_codon + '-to-' + new_codon
															
															if silent_mutation_relative_to_edit == 'overlap': 
																if 'PAM_disrupted' in pe_annotate: 
																	pe_annotate = 'silent_mutation_and_PAM_disrupted'
																else:
																	pe_annotate = 'silent_mutation_and_PAM_intact'
															else:
																pe_annotate = 'silent_mutation_and_PAM_disrupted'
																
															silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
															silent_mutation_relative_to_edit = 'upstream'
															break
												
												else: # If methionine or tryptophan codons, try to introduce silent mutation 3' of the PAM
													codon_start_idx = 7
													codon_end_idx = 10
													original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
													aa_identity = codon_dict[original_codon][1]

													if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
														for codon_substitute in aa2codon[aa_identity]:
															if codon_substitute[0] != original_codon:

																new_codon = codon_substitute[0].lower()
																pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + rtt_length])
																pe_silent_mutation = original_codon + '-to-' + new_codon

																if silent_mutation_relative_to_edit == 'overlap': 
																	if 'PAM_disrupted' in pe_annotate: 
																		pe_annotate = 'PAM_disrupted_and_silent_mutation'
																	else:
																		pe_annotate = 'PAM_intact_and_silent_mutation'
																else:
																	pe_annotate = 'PAM_disrupted_and_silent_mutation'
																
																silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
																silent_mutation_relative_to_edit = 'downstream'
																break
								
							# Check to see if pegRNA extension is within input sequence
							if len(pegRNA_ext) == (pbs_length + rtt_length):

								# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
								if pegid not in pe_design[target_name]:

									# First list is for peg extension, second list is for nicking guide
									pe_design[target_name][pegid] = [[],[]]

								# Store pegRNA design
								if silent_mutation_edit_sequence == '':
									if pe_pam_ref_silent_mutation == '':
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])

									else:
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref_silent_mutation, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])
								else:
									if pe_pam_ref_silent_mutation == '':
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit])

									else:
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref_silent_mutation, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit])

				# Create ngRNAs targeting (-) strand for (+) pegRNAs
				if pegid in pe_design[target_name]:
					for ng_minus in target_design[target_name]['ngRNA']['-']:
						ng_nick_ref_idx, ng_edit_start_idx, ng_edit_end_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_minus
						nick_distance = ng_nick_ref_idx - pe_nick_ref_idx

						if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]') and (len(silent_mutation_edit_sequence) > 0):
							ng_spacer_sequence_edit = silent_mutation_edit_sequence[ng_edit_start_idx:ng_edit_end_idx]

							mutation_indices = [i for i, a in enumerate(ng_spacer_sequence_edit) if a.islower()]
							if len(mutation_indices) > 0:
								if len([1 for x in mutation_indices if x < 10]) > 0:
									ng_annotate = 'PE3b-seed'

								else:
									ng_annotate = 'PE3b-nonseed'
							else:
								ng_annotate = 'PE3'

						if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):
							pe_design[target_name][pegid][1].append([ng_nick_ref_idx, reverse_complement(ng_spacer_sequence_edit), reverse_complement(ng_pam_edit), ng_annotate, '-', nick_distance])

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
							if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]'):
								
								if pe_annotate_constant == 'PAM_intact':
									
									if nick_aa_index == 0:
										codon_start_idx = -6
										codon_end_idx = -3
										original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
										aa_identity = codon_dict[original_codon][1]

										for codon_substitute in aa2codon[aa_identity]:
											if not re.search(reverse_complement(pam_search), codon_substitute[0]):

												new_codon = codon_substitute[0].lower()
												pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
												pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(new_codon)
												pe_annotate = 'PAM_disrupted_with_silent_mutation'
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
												if nick2lastedit_length >= abs(codon_start_idx): silent_mutation_relative_to_edit = 'downstream'
												elif nick2edit_length <= abs(codon_end_idx): silent_mutation_relative_to_edit = 'upstream'
												else: silent_mutation_relative_to_edit = 'overlap'
												break
										
										if pe_annotate != 'PAM_disrupted_with_silent_mutation': # Generate silent mutation that does not disrupt PAM sequence
									
											for codon_substitute in aa2codon[aa_identity]:
												if codon_substitute[0] != original_codon:

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
													pe_silent_mutation = original_codon + '-to-' + new_codon
													pe_annotate = 'PAM_intact_with_silent_mutation'
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
													if nick2lastedit_length >= abs(codon_start_idx): silent_mutation_relative_to_edit = 'downstream'
													elif nick2edit_length <= abs(codon_end_idx): silent_mutation_relative_to_edit = 'upstream'
													else: silent_mutation_relative_to_edit = 'overlap'
													break

									elif nick_aa_index == 1:
										codon_start_idx = -7
										codon_end_idx = -4
										original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
										aa_identity = codon_dict[original_codon][1]

										for codon_substitute in aa2codon[aa_identity]:

											pam_slice = codon_substitute[0][1:] + edit_sequence[pe_nick_edit_idx - 4:pe_nick_edit_idx - 3].upper()
											if not re.search(reverse_complement(pam_search), pam_slice):

												new_codon = codon_substitute[0].lower()
												pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
												pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(pam_slice).lower()
												pe_annotate = 'PAM_disrupted_with_silent_mutation'
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
												if nick2lastedit_length >= abs(codon_start_idx): silent_mutation_relative_to_edit = 'downstream'
												elif nick2edit_length <= abs(codon_end_idx): silent_mutation_relative_to_edit = 'upstream'
												else: silent_mutation_relative_to_edit = 'overlap'
												break

										if pe_annotate != 'PAM_disrupted_with_silent_mutation': # Generate silent mutation that does not disrupt PAM sequence
									
											for codon_substitute in aa2codon[aa_identity]:
												if codon_substitute[0] != original_codon:

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
													pe_silent_mutation = original_codon + '-to-' + new_codon
													pe_annotate = 'PAM_intact_with_silent_mutation'
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
													if nick2lastedit_length >= abs(codon_start_idx): silent_mutation_relative_to_edit = 'downstream'
													elif nick2edit_length <= abs(codon_end_idx): silent_mutation_relative_to_edit = 'upstream'
													else: silent_mutation_relative_to_edit = 'overlap'
													break

									elif nick_aa_index == 2:
										codon_start_idx_1 = -8
										codon_end_idx_1 = -5
										codon_start_idx_2 = -5
										codon_end_idx_2 = -2

										original_codon_1 = edit_sequence[pe_nick_edit_idx + codon_start_idx_1:pe_nick_edit_idx + codon_end_idx_1].upper()
										original_codon_2 = edit_sequence[pe_nick_edit_idx + codon_start_idx_2:pe_nick_edit_idx + codon_end_idx_2].upper()

										aa_identity_1 = codon_dict[original_codon_1][1]
										aa_identity_2 = codon_dict[original_codon_2][1]

										for codon_substitute in aa2codon[aa_identity_1]:

											new_codons = codon_substitute[0] + original_codon_2
											pam_slice = new_codons[2:5]

											if not re.search(reverse_complement(pam_search), pam_slice):

												new_codon = codon_substitute[0].lower()
												pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + pbs_length]
												pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(pam_slice).lower()
												pe_annotate = 'PAM_disrupted_with_silent_mutation'
												silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:]
												if nick2lastedit_length >= abs(codon_start_idx_1): silent_mutation_relative_to_edit = 'downstream'
												elif nick2edit_length <= abs(codon_end_idx_1): silent_mutation_relative_to_edit = 'upstream'
												else: silent_mutation_relative_to_edit = 'overlap'
												break

										if pe_annotate != 'PAM_disrupted_with_silent_mutation':

											for codon_substitute in aa2codon[aa_identity_2]:

												new_codons = original_codon_1 + codon_substitute[0]
												pam_slice = new_codons[2:5]

												if not re.search(reverse_complement(pam_search), pam_slice):

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + pbs_length]
													pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(pam_slice).lower()
													pe_annotate = 'PAM_disrupted_with_silent_mutation'
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:]
													if nick2lastedit_length >= abs(codon_start_idx_2): silent_mutation_relative_to_edit = 'downstream'
													elif nick2edit_length <= abs(codon_end_idx_2): silent_mutation_relative_to_edit = 'upstream'
													else: silent_mutation_relative_to_edit = 'overlap'
													break
										
										if 'silent_mutation' not in pe_annotate: # Generate silent mutation that does not disrupt PAM sequence
											
											for codon_substitute in aa2codon[aa_identity_1]: # Try to introduce silent mutation in codon 1 first
												if codon_substitute[0] != original_codon_1:

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:pe_nick_edit_idx + pbs_length]
													pe_silent_mutation = original_codon_1 + '-to-' + new_codon
													pe_annotate = 'PAM_intact_with_silent_mutation'
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_1] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_1:]
													if nick2lastedit_length >= abs(codon_start_idx_1): silent_mutation_relative_to_edit = 'downstream'
													elif nick2edit_length <= abs(codon_end_idx_1): silent_mutation_relative_to_edit = 'upstream'
													else: silent_mutation_relative_to_edit = 'overlap'
													break
											
										if 'silent_mutation' not in pe_annotate: # Try to introduce silent mutation in codon 2 if codon 1 did not work

											for codon_substitute in aa2codon[aa_identity_2]:
												if codon_substitute[0] != original_codon_2:

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:pe_nick_edit_idx + pbs_length]
													pe_silent_mutation = original_codon_2 + '-to-' + new_codon
													pe_annotate = 'PAM_intact_with_silent_mutation'
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx_2] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx_2:]
													if nick2lastedit_length >= abs(codon_start_idx_2): silent_mutation_relative_to_edit = 'downstream'
													elif nick2edit_length <= abs(codon_end_idx_2): silent_mutation_relative_to_edit = 'upstream'
													else: silent_mutation_relative_to_edit = 'overlap'
													break

								if pe_annotate_constant == 'PAM_disrupted' or silent_mutation_relative_to_edit == 'overlap': # Generate silent mutation 5' or 3' of the PAM
							
									if nick_aa_index == 0:
										codon_start_idx = -3 # Try to introduce silent mutation 3' of the PAM first
										codon_end_idx = 0
										original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
										aa_identity = codon_dict[original_codon][1]

										if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
											for codon_substitute in aa2codon[aa_identity]:
												if codon_substitute[0] != original_codon:

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
													pe_silent_mutation = original_codon + '-to-' + new_codon
													
													if silent_mutation_relative_to_edit == 'overlap': 
														if 'PAM_disrupted' in pe_annotate: 
															pe_annotate = 'PAM_disrupted_and_silent_mutation'
														else:
															pe_annotate = 'PAM_intact_and_silent_mutation'
													else:
														pe_annotate = 'PAM_disrupted_and_silent_mutation'
													
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
													silent_mutation_relative_to_edit = 'downstream'
													break

										else: # If methionine or tryptophan codons, try to introduce silent mutation 5' of the PAM
											codon_start_idx = -9
											codon_end_idx = -6
											original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
											aa_identity = codon_dict[original_codon][1]

											if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
												for codon_substitute in aa2codon[aa_identity]:
													if codon_substitute[0] != original_codon:
														
														new_codon = codon_substitute[0].lower()
														pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
														pe_silent_mutation = original_codon + '-to-' + new_codon
												
														if silent_mutation_relative_to_edit == 'overlap': 
															if 'PAM_disrupted' in pe_annotate: 
																pe_annotate = 'silent_mutation_and_PAM_disrupted'
															else:
																pe_annotate = 'silent_mutation_and_PAM_intact'
														else:
															pe_annotate = 'silent_mutation_and_PAM_disrupted'
														
														silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
														silent_mutation_relative_to_edit = 'upstream'
														break
									
									elif nick_aa_index == 1:
										codon_start_idx = -4 # Try to introduce silent mutation 3' of the PAM first
										codon_end_idx = -1
										original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
										aa_identity = codon_dict[original_codon][1]

										if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
											for codon_substitute in aa2codon[aa_identity]:
												if codon_substitute[0] != original_codon:

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
													pe_silent_mutation = original_codon + '-to-' + new_codon

													if silent_mutation_relative_to_edit == 'overlap': 
														if 'PAM_disrupted' in pe_annotate: 
															pe_annotate = 'PAM_disrupted_and_silent_mutation'
														else:
															pe_annotate = 'PAM_intact_and_silent_mutation'
													else:
														pe_annotate = 'PAM_disrupted_and_silent_mutation'

													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
													silent_mutation_relative_to_edit = 'downstream'
													break

										else: # If methionine or tryptophan codons, try to introduce silent mutation 5' of the PAM
											codon_start_idx = -10
											codon_end_idx = -7
											original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
											aa_identity = codon_dict[original_codon][1]

											if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
												for codon_substitute in aa2codon[aa_identity]:
													if codon_substitute[0] != original_codon:
														
														new_codon = codon_substitute[0].lower()
														pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
														pe_silent_mutation = original_codon + '-to-' + new_codon
														
														if silent_mutation_relative_to_edit == 'overlap': 
															if 'PAM_disrupted' in pe_annotate: 
																pe_annotate = 'silent_mutation_and_PAM_disrupted'
															else:
																pe_annotate = 'silent_mutation_and_PAM_intact'
														else:
															pe_annotate = 'silent_mutation_and_PAM_disrupted'
														
														silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
														silent_mutation_relative_to_edit = 'upstream'
														break
									
									elif nick_aa_index == 2:
										codon_start_idx = -2 # Try to introduce silent mutation 3' of the PAM first
										codon_end_idx = 1
										original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
										aa_identity = codon_dict[original_codon][1]

										if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
											for codon_substitute in aa2codon[aa_identity]:
												if (codon_substitute[0][2] == original_codon[2]) & (codon_substitute[0] != original_codon): # Keep the third base the same to preserve PBS sequence

													new_codon = codon_substitute[0].lower()
													pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
													pe_silent_mutation = original_codon + '-to-' + new_codon
											
													if silent_mutation_relative_to_edit == 'overlap': 
														if 'PAM_disrupted' in pe_annotate: 
															pe_annotate = 'PAM_disrupted_and_silent_mutation'
														else:
															pe_annotate = 'PAM_intact_and_silent_mutation'
													else:
														pe_annotate = 'PAM_disrupted_and_silent_mutation'
											
													silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
													silent_mutation_relative_to_edit = 'downstream'
													break

										else: # If methionine or tryptophan codons, try to introduce silent mutation 5' of the PAM
											codon_start_idx = -11
											codon_end_idx = -8
											original_codon = edit_sequence[pe_nick_edit_idx + codon_start_idx:pe_nick_edit_idx + codon_end_idx].upper()
											aa_identity = codon_dict[original_codon][1]

											if (original_codon != 'ATG') & (original_codon != 'TGG'): # Can not introduce silent mutations for methionine and tryptophan codons
												for codon_substitute in aa2codon[aa_identity]:
													if codon_substitute[0] != original_codon:
														
														new_codon = codon_substitute[0].lower()
														pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:pe_nick_edit_idx + pbs_length]
														pe_silent_mutation = original_codon + '-to-' + new_codon
												
														if silent_mutation_relative_to_edit == 'overlap': 
															if 'PAM_disrupted' in pe_annotate: 
																pe_annotate = 'silent_mutation_and_PAM_disrupted'
															else:
																pe_annotate = 'silent_mutation_and_PAM_intact'
														else:
															pe_annotate = 'silent_mutation_and_PAM_disrupted'
														
														silent_mutation_edit_sequence = edit_sequence[:pe_nick_edit_idx + codon_start_idx] + new_codon + edit_sequence[pe_nick_edit_idx + codon_end_idx:]
														silent_mutation_relative_to_edit = 'upstream'
														break

							# Check to see if pegRNA extension is within input sequence
							if len(pegRNA_ext) == (pbs_length + rtt_length):

								# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
								if pegid not in pe_design[target_name]:

									# First list is for peg extension, second list is for nicking guide
									pe_design[target_name][pegid] = [[],[]]

								# Store pegRNA design
								if silent_mutation_edit_sequence == '':
									if pe_pam_ref_silent_mutation == '':
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), reverse_complement(pe_pam_ref), pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])
									
									else:
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), pe_pam_ref_silent_mutation, pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])
								
								else:
									if pe_pam_ref_silent_mutation == '':
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), reverse_complement(pe_pam_ref), pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit])
									
									else:
										pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), pe_pam_ref_silent_mutation, pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, silent_mutation_edit_sequence, silent_mutation_relative_to_edit])

				# Create ngRNAs targeting (+) strand for (-) pegRNAs
				if pegid in pe_design[target_name]:
					for ng_plus in target_design[target_name]['ngRNA']['+']:
						ng_nick_ref_idx, ng_edit_start_idx, ng_edit_end_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_plus
						nick_distance = ng_nick_ref_idx - pe_nick_ref_idx

						if silent_mutation and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]') and (len(silent_mutation_edit_sequence) > 0):
							ng_spacer_sequence_edit = silent_mutation_edit_sequence[ng_edit_start_idx:ng_edit_end_idx]

							mutation_indices = [i for i, a in enumerate(ng_spacer_sequence_edit) if a.islower()]
							if len(mutation_indices) > 0:
								if len([1 for x in mutation_indices if x >= 10]) > 0:
									ng_annotate = 'PE3b-seed'

								else:
									ng_annotate = 'PE3b-nonseed'
							else:
								ng_annotate = 'PE3'

						if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):
							pe_design[target_name][pegid][1].append([ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, '+', nick_distance])

		if counter%1000 == 0:
			logger.info('Completed pegRNA and ngRNA search for %s out of %s sites ...' % (counter, total_regions))
		counter += 1

logger.info('Completed pegRNA and ngRNA search for %s out of %s sites ...' % (counter - 1, total_regions))

# Output pegRNAs
pegRNAs_summary_f = '%s_PrimeDesign.csv' % str(time.strftime("%Y%m%d_%I.%M.%S", time.localtime()))
logger.info('Writing pegRNA and ngRNA designs into output file %s ...' % pegRNAs_summary_f)

counter = 1
with open(out_dir + '/%s' % pegRNAs_summary_f, 'w') as f:

	f.write(','.join(map(str, ['Target_name', 'Target_sequence', 'pegRNA_number', 'gRNA_type', 'Spacer_sequence', 'Spacer_GC_content', 'PAM_sequence', 'Extension_sequence', 'Strand', 'Annotation', 'pegRNA-to-edit_distance', 'Nick_index', 'ngRNA-to-pegRNA_distance', 'PBS_length', 'PBS_GC_content', 'RTT_length', 'RTT_GC_content', 'First_extension_nucleotide', 'Spacer_sequence_order_TOP', 'Spacer_sequence_order_BOTTOM', 'pegRNA_extension_sequence_order_TOP', 'pegRNA_extension_sequence_order_BOTTOM', 'Edit_type', 'Reference_sequence', 'Edit_sequence','Silent_mutation_relative_to_edit'])) + '\n')
	for target_name in pe_design:

		if genome_wide_design or saturation_mutagenesis:

			if filter_c1_extension:
				
				peg_count = 0
				for pegid in list(pe_design[target_name].keys()):

					ng_continue = True

					# Write pegRNAs
					for pegRNA_entry in pe_design[target_name][pegid][0]:
						pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, pe_strand, pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit = pegRNA_entry

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
										f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')
										peg_count += 1

									else:
										ng_continue = False

								else:

									f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')
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
									f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')
									peg_count += 1
									
								else:
									ng_continue = False

							else:

								f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')
								peg_count += 1

						# Write ngRNAs
						if ng_continue:
							for ngRNA_entry in pe_design[target_name][pegid][1][:number_of_ngrnas]:
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
										f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')

								else:
									f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')

						counter += 1

						if peg_count == number_of_pegrnas:
							break

			else:

				for pegid in list(pe_design[target_name].keys())[:number_of_pegrnas]:

					# Write pegRNAs
					for pegRNA_entry in pe_design[target_name][pegid][0]:
						pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, pe_strand, pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit = pegRNA_entry

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
								f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')

						else:

							f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')

						# Sort ngRNAs

						# Write ngRNAs
						for ngRNA_entry in pe_design[target_name][pegid][1][:number_of_ngrnas]:
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
									f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')

							else:
								f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')

						counter += 1

		else:
		
			for pegid in pe_design[target_name]:

				ng_continue = True

				# Write pegRNAs
				for pegRNA_entry in pe_design[target_name][pegid][0]:
					pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, pe_strand, pbs_length, rtt_length, pegRNA_ext, nick2lastedit_length, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit = pegRNA_entry

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
							f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')

						else:
							ng_continue = False

					else:

						f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, spacer_oligo_top, spacer_oligo_bottom, pegext_oligo_top, pegext_oligo_bottom, edit_type, reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')

				# Write ngRNAs
				if ng_continue:
					for ngRNA_entry in pe_design[target_name][pegid][1][:number_of_ngrnas]:
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
								f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')

						else:

							f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', spacer_oligo_top, spacer_oligo_bottom, '', '', '', reference_sequence, edit_sequence, silent_mutation_relative_to_edit])) + '\n')

				counter += 1