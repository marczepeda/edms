#!/usr/bin/env python3
"""Design edit sequences for deep mutational scanning.

DMSDesign is a lightweight PrimeDesign-style script for generating edited DNA
sequences, including saturation mutagenesis libraries. It does not design
pegRNAs, spacers, PBSs, RTTs, PAM edits, or nicking guides.
"""

import argparse
import csv
import difflib
import logging
import os
import re
import sys
import time
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='''----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Software for the design of edit sequences for deep mutational scanning.
Input sequence syntax follows PrimeDesign: substitutions (REF/EDIT), insertions (+EDIT), deletions (-REF).
For saturation mutagenesis, place one region in parentheses and choose -sat_mut.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------''', formatter_class=RawTextHelpFormatter)

parser.add_argument('-f', '--file', required=True, type=str, help="""Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (Required)

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

""")
parser.add_argument('-sat_mut', '--saturation_mutagenesis', default=False, choices=['aa', 'aa_subs', 'aa_ins', 'aa_dels', 'aa_silent', 'base'], type=str,
                    help='Saturation mutagenesis mode. aa = aa_subs + aa_ins + aa_dels. Default: False')
parser.add_argument('-silent_mut', '--silent_mutation', default=0, type=int,
                    help='Number of additional silent codon changes to introduce into each edit sequence. Must be >= 0. Default: 0')
parser.add_argument('-silent_mut_mode', '--silent_mutation_mode', default='barcode',
                    choices=['close', 'upstream', 'downstream', 'distribute', 'barcode'], type=str,
                    help='Silent mutation placement mode. close = nearest intended edit; upstream = 5-prime of edit; downstream = 3-prime of edit; distribute = spread across allowed region; barcode = deterministic variant-specific synonymous pattern. Default: barcode')
parser.add_argument('-out', '--out_dir', default='./DMSDesign/DATETIMESTAMP_DMSDesign', type=str, help='Output directory. Default: ./DMSDesign/DATETIMESTAMP_DMSDesign')

args = parser.parse_args()
file_in = args.file
saturation_mutagenesis = args.saturation_mutagenesis
silent_mutation = args.silent_mutation
silent_mutation_mode = args.silent_mutation_mode

if silent_mutation < 0:
    raise ValueError('--silent_mutation must be an integer greater than or equal to 0')

out_dir = args.out_dir
if out_dir == './DMSDesign/DATETIMESTAMP_DMSDesign':
    out_dir = './DMSDesign/%s_DMSDesign' % str(time.strftime('%y%m%d_%H.%M.%S', time.localtime()))
os.makedirs(out_dir, exist_ok=True)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh = logging.FileHandler(os.path.join(out_dir, 'DMSDesign.log'))
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

codon_dict = {
    'GGG':['Gly','G',0.25],'GGA':['Gly','G',0.25],'GGT':['Gly','G',0.16],'GGC':['Gly','G',0.34],
    'GAG':['Glu','E',0.58],'GAA':['Glu','E',0.42],'GAT':['Asp','D',0.46],'GAC':['Asp','D',0.54],
    'GTG':['Val','V',0.47],'GTA':['Val','V',0.11],'GTT':['Val','V',0.18],'GTC':['Val','V',0.24],
    'GCG':['Ala','A',0.11],'GCA':['Ala','A',0.23],'GCT':['Ala','A',0.26],'GCC':['Ala','A',0.4],
    'AGG':['Arg','R',0.2],'AGA':['Arg','R',0.2],'AGT':['Ser','S',0.15],'AGC':['Ser','S',0.24],
    'AAG':['Lys','K',0.58],'AAA':['Lys','K',0.42],'AAT':['Asn','N',0.46],'AAC':['Asn','N',0.54],
    'ATG':['Met','M',1],'ATA':['Ile','I',0.16],'ATT':['Ile','I',0.36],'ATC':['Ile','I',0.48],
    'ACG':['Thr','T',0.12],'ACA':['Thr','T',0.28],'ACT':['Thr','T',0.24],'ACC':['Thr','T',0.36],
    'TGG':['Trp','W',1],'TGA':['End','X',0.52],'TGT':['Cys','C',0.45],'TGC':['Cys','C',0.55],
    'TAG':['End','X',0.2],'TAA':['End','X',0.28],'TAT':['Tyr','Y',0.43],'TAC':['Tyr','Y',0.57],
    'TTG':['Leu','L',0.13],'TTA':['Leu','L',0.07],'TTT':['Phe','F',0.45],'TTC':['Phe','F',0.55],
    'TCG':['Ser','S',0.06],'TCA':['Ser','S',0.15],'TCT':['Ser','S',0.18],'TCC':['Ser','S',0.22],
    'CGG':['Arg','R',0.21],'CGA':['Arg','R',0.11],'CGT':['Arg','R',0.08],'CGC':['Arg','R',0.19],
    'CAG':['Gln','Q',0.75],'CAA':['Gln','Q',0.25],'CAT':['His','H',0.41],'CAC':['His','H',0.59],
    'CTG':['Leu','L',0.41],'CTA':['Leu','L',0.07],'CTT':['Leu','L',0.13],'CTC':['Leu','L',0.2],
    'CCG':['Pro','P',0.11],'CCA':['Pro','P',0.27],'CCT':['Pro','P',0.28],'CCC':['Pro','P',0.33],
}

aa2codon = {}
for codon, values in codon_dict.items():
    aa2codon.setdefault(values[1], []).append([codon, values[2]])
for aa in aa2codon:
    aa2codon[aa] = sorted(aa2codon[aa], key=lambda x: x[1], reverse=True)


def codon_nt_difference(codon_ref, codon_edit):
    return sum(1 for ref_base, edit_base in zip(codon_ref.upper(), codon_edit.upper()) if ref_base != edit_base)


def strip_saturation_parentheses(target_sequence):
    start = target_sequence.find('(')
    stop = target_sequence.find(')')
    sequence_left = target_sequence[:start]
    sequence_to_edit = target_sequence[start + 1:stop]
    sequence_right = target_sequence[stop + 1:]
    return sequence_left, sequence_to_edit, sequence_right


def validate_saturation_input(target_sequence):
    format_check = ''.join([i for i in target_sequence if i in ['(', ')', '/', '+', '-']])
    if len(target_sequence) != sum(1 for x in target_sequence.upper() if x in ['A', 'T', 'C', 'G', '(', ')']):
        logger.error('Input sequence %s contains invalid characters for saturation mutagenesis.', target_sequence)
        sys.exit(1)
    if format_check.count('(') != 1 or format_check.count(')') != 1:
        logger.error('Input sequence %s must have exactly one complete parenthesized saturation region.', target_sequence)
        sys.exit(1)


def process_sequence(input_sequence):
    format_check = ''.join([i for i in input_sequence if i in ['(', ')', '/', '+', '-']])
    if len(input_sequence) != sum(1 for x in input_sequence.upper() if x in ['A', 'T', 'C', 'G', '(', ')', '+', '-', '/']):
        logger.error('Input sequence %s contains invalid characters.', input_sequence)
        sys.exit(1)
    if not (format_check.count('(') == format_check.count(')') and format_check.count('(') > 0):
        logger.error('Input sequence %s does not have full sets of parentheses.', input_sequence)
        sys.exit(1)
    if '((' in format_check or '()' in format_check:
        logger.error('Input sequence %s has nested or empty parentheses.', input_sequence)
        sys.exit(1)
    if sum(1 if x in format_check else 0 for x in ['++','--','//','+-','+/','-+','-/','/+','/-','/(','+(','-(',')/',')+',')-']) != 0:
        logger.error('Input sequence %s has more than one edit annotation per parentheses set.', input_sequence)
        sys.exit(1)

    editformat2sequence = {}
    editnumber2sequence = {}
    edit_idxs = [[m.start(), m.end()] for m in re.finditer(r'\(.*?\)', input_sequence)]
    edit_counter = 1
    for edit_idx in edit_idxs:
        edit = input_sequence[edit_idx[0]:edit_idx[1]]
        if '/' in edit:
            ref = edit.split('/')[0].replace('(', '')
            alt = edit.split('/')[1].replace(')', '').lower()
        elif '+' in edit:
            ref = ''
            alt = edit.split('+')[1].replace(')', '').lower()
        elif '-' in edit:
            ref = edit.split('-')[1].replace(')', '')
            alt = ''
        editformat2sequence[edit] = [ref, alt, edit_counter]
        editnumber2sequence[edit_counter] = [ref, alt]
        edit_counter += 1

    edit_start = min(i.start() for i in re.finditer(r'\(', input_sequence))
    edit_stop = max(i.start() for i in re.finditer(r'\)', input_sequence))
    edit_span_sequence_w_ref = input_sequence[edit_start:edit_stop + 1]
    edit_span_sequence_w_edit = input_sequence[edit_start:edit_stop + 1]
    reference_sequence = input_sequence
    edit_sequence = input_sequence
    editnumber_sequence = input_sequence

    for edit, values in editformat2sequence.items():
        edit_span_sequence_w_ref = edit_span_sequence_w_ref.replace(edit, values[0])
        edit_span_sequence_w_edit = edit_span_sequence_w_edit.replace(edit, values[1])
        reference_sequence = reference_sequence.replace(edit, values[0])
        edit_sequence = edit_sequence.replace(edit, values[1])
        editnumber_sequence = editnumber_sequence.replace(edit, str(values[2]))

    edit_start_in_ref = re.search(r'\(', input_sequence).start()
    edit_stop_in_ref_rev = re.search(r'\)', input_sequence[::-1]).start()
    return {
        'editformat2sequence': editformat2sequence,
        'editnumber2sequence': editnumber2sequence,
        'reference_sequence': reference_sequence,
        'edit_sequence': edit_sequence,
        'editnumber_sequence': editnumber_sequence,
        'edit_span_length_w_ref': len(edit_span_sequence_w_ref),
        'edit_span_length_w_edit': len(edit_span_sequence_w_edit),
        'edit_start_in_ref': edit_start_in_ref,
        'edit_stop_in_ref_rev': edit_stop_in_ref_rev,
    }


def edit_type_from_target(target_sequence):
    edit_type = ''
    if '/' in target_sequence:
        edit_type += '& substitution'
    if '+' in target_sequence:
        edit_type += '& insertion'
    if '-' in target_sequence:
        edit_type += '& deletion'
    return edit_type[2:]


def changed_positions(reference_sequence, edit_sequence):
    positions = [m.start() for m in re.finditer('[acgt]', edit_sequence)]
    if positions:
        return positions
    # fallback for deletions, where edited bases may be absent
    sm = difflib.SequenceMatcher(None, reference_sequence.upper(), edit_sequence.upper())
    out = []
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag != 'equal':
            out.extend(range(j1, max(j2, j1 + 1)))
    return out


def candidate_silent_codons(edit_sequence, reference_sequence, mutation_positions, allowed_range=None, mode='close'):
    seq = edit_sequence.upper()
    candidates = []
    if not mutation_positions:
        mutation_center = len(seq) / 2
        mutation_min = mutation_center
        mutation_max = mutation_center
    else:
        mutation_center = sum(mutation_positions) / len(mutation_positions)
        mutation_min = min(mutation_positions)
        mutation_max = max(mutation_positions)

    for codon_start in range(0, len(seq) - 2, 3):
        codon_end = codon_start + 3
        codon_center = codon_start + 1

        if allowed_range is not None:
            if codon_start < allowed_range[0] or codon_end > allowed_range[1]:
                continue

        if mode == 'upstream' and codon_end > mutation_min:
            continue
        if mode == 'downstream' and codon_start <= mutation_max:
            continue

        codon_ref = seq[codon_start:codon_end]
        if codon_ref not in codon_dict:
            continue
        if any(codon_start <= pos < codon_end for pos in mutation_positions):
            continue
        aa = codon_dict[codon_ref][1]
        synonymous = [x for x in aa2codon[aa] if x[0] != codon_ref]
        if not synonymous:
            continue
        synonymous = sorted(synonymous, key=lambda x: (-codon_nt_difference(codon_ref, x[0]), -x[1], x[0]))
        distance = abs(codon_center - mutation_center)
        candidates.append((distance, codon_start, codon_end, codon_ref, synonymous))

    if mode == 'upstream':
        return sorted(candidates, key=lambda x: (x[0], -x[1]))
    if mode == 'downstream':
        return sorted(candidates, key=lambda x: (x[0], x[1]))
    return sorted(candidates, key=lambda x: (x[0], x[1]))



def choose_distributed_candidates(candidates, n_silent, allowed_range, edit_sequence):
    if n_silent <= 0 or not candidates:
        return []
    if allowed_range is None:
        start = 0
        stop = len(edit_sequence)
    else:
        start, stop = allowed_range
    if n_silent == 1:
        targets = [(start + stop - 1) / 2]
    else:
        targets = [start + ((stop - start - 1) * i / (n_silent - 1)) for i in range(n_silent)]

    selected = []
    used_starts = set()
    for target in targets:
        available = [c for c in candidates if c[1] not in used_starts]
        if not available:
            break
        chosen = min(available, key=lambda c: (abs((c[1] + 1) - target), c[0], c[1]))
        selected.append(chosen)
        used_starts.add(chosen[1])
    return selected


def add_silent_mutations(edit_sequence, reference_sequence, n_silent, allowed_range=None, mode='close', barcode_seed=0):
    if n_silent == 0:
        return edit_sequence, '', '', 0
    if mode not in ['close', 'upstream', 'downstream', 'distribute']:
        raise ValueError("mode must be one of: 'close', 'upstream', 'downstream', 'distribute'; 'barcode' mode should use apply_silent_barcode() instead")

    edited = edit_sequence
    mutation_positions = changed_positions(reference_sequence, edit_sequence)
    records = []

    if mode == 'distribute':
        candidates = candidate_silent_codons(edited, reference_sequence, mutation_positions, allowed_range=allowed_range, mode='close')
        selected_candidates = choose_distributed_candidates(candidates, n_silent, allowed_range, edited)
        for _, codon_start, codon_end, codon_ref, synonymous in selected_candidates:
            codon_edit = synonymous[0][0].lower()
            edited = edited[:codon_start] + codon_edit + edited[codon_end:]
            records.append({
                'codon_start': codon_start,
                'codon_end': codon_end,
                'codon_ref': codon_ref,
                'codon_edit': codon_edit,
                'aa': codon_dict[codon_ref][1],
            })
    else:
        for silent_idx in range(n_silent):
            candidates = candidate_silent_codons(edited, reference_sequence, mutation_positions, allowed_range=allowed_range, mode=mode)
            candidates = [c for c in candidates if c[1] not in [r['codon_start'] for r in records]]
            if not candidates:
                break
            _, codon_start, codon_end, codon_ref, synonymous = candidates[0]
            codon_edit = synonymous[0][0].lower()
            edited = edited[:codon_start] + codon_edit + edited[codon_end:]
            records.append({
                'codon_start': codon_start,
                'codon_end': codon_end,
                'codon_ref': codon_ref,
                'codon_edit': codon_edit,
                'aa': codon_dict[codon_ref][1],
            })

    annotation = ';'.join('%s%d%s>%s' % (r['aa'], int(r['codon_start'] / 3) + 1, r['codon_ref'], r['codon_edit']) for r in records)
    positions = ';'.join('%d-%d' % (r['codon_start'], r['codon_end'] - 1) for r in records)
    return edited, annotation, positions, len(records)

def barcode_candidate_pool(edit_sequence, reference_sequence, n_silent, allowed_range=None):
    """
    Build the silent-codon barcode channel.

    Uses the same positional constraint as close mode:
    eligible codons are ranked by distance from the intended edit.
    """
    if n_silent <= 0:
        return []

    mutation_positions = changed_positions(reference_sequence, edit_sequence)

    candidates = candidate_silent_codons(
        edit_sequence=edit_sequence,
        reference_sequence=reference_sequence,
        mutation_positions=mutation_positions,
        allowed_range=allowed_range,
        mode='close',
    )

    pool = []
    for _distance, codon_start, codon_end, codon_ref, synonymous in candidates[:n_silent]:
        choices = sorted(
            synonymous,
            key=lambda x: (-codon_nt_difference(codon_ref, x[0]), -x[1], x[0])
        )

        pool.append({
            'codon_start': codon_start,
            'codon_end': codon_end,
            'codon_ref': codon_ref,
            'aa': codon_dict[codon_ref][1],
            'choices': [choice[0] for choice in choices],
        })

    return pool


def barcode_capacity(pool):
    """
    Number of unique barcode patterns possible from the selected silent codons.
    """
    if not pool:
        return 0

    capacity = 1
    for candidate in pool:
        capacity *= len(candidate['choices'])

    return capacity


def cyclic_barcode_indices(barcode_index, pool):
    """
    Convert a barcode index into one synonymous-codon choice per barcode codon.

    This is mixed-radix counting, so choices cycle evenly:
    position 1 changes fastest, then position 2, etc.
    """
    indices = []

    for candidate in pool:
        n_choices = len(candidate['choices'])
        indices.append(barcode_index % n_choices)
        barcode_index //= n_choices

    return indices

def apply_silent_barcode(edit_sequence, reference_sequence, n_silent, barcode_id, allowed_range=None):
    pool = barcode_candidate_pool(
        edit_sequence,
        reference_sequence,
        n_silent,
        allowed_range=allowed_range,
    )

    capacity = barcode_capacity(pool)
    if capacity == 0:
        return edit_sequence, '', '', 0, '', '', ''

    barcode_index = barcode_id % capacity
    reuse_count = barcode_id // capacity

    edited = edit_sequence
    records = []

    choice_indices = cyclic_barcode_indices(barcode_index, pool)

    for candidate, choice_index in zip(pool, choice_indices):
        codon_edit = candidate["choices"][choice_index].lower()
        codon_start = candidate["codon_start"]
        codon_end = candidate["codon_end"]

        edited = edited[:codon_start] + codon_edit + edited[codon_end:]

        records.append({
            "codon_start": codon_start,
            "codon_end": codon_end,
            "codon_ref": candidate["codon_ref"],
            "codon_edit": codon_edit,
            "aa": candidate["aa"],
        })

    annotation = ";".join(
        "%s%d%s>%s" % (
            r["aa"],
            int(r["codon_start"] / 3) + 1,
            r["codon_ref"],
            r["codon_edit"],
        )
        for r in records
    )

    positions = ";".join(
        "%d-%d" % (r["codon_start"], r["codon_end"] - 1)
        for r in records
    )

    pattern = ";".join(
        "%d:%s>%s" % (
            int(r["codon_start"] / 3) + 1,
            r["codon_ref"],
            r["codon_edit"],
        )
        for r in records
    )

    return edited, annotation, positions, len(records), pattern, barcode_index, reuse_count

def sorted_codons_by_difference(codon_ref, codon_edit_entry_list):
    return sorted(codon_edit_entry_list, key=lambda x: (-codon_nt_difference(codon_ref, x[0]), -x[1], x[0]))


def sorted_codons_by_usage(codon_edit_entry_list):
    return sorted(codon_edit_entry_list, key=lambda x: (-x[1], x[0]))


def saturating_mutagenesis_input_sequences(target_name, target_sequence, sm_type):
    validate_saturation_input(target_sequence)
    sequence_left, sequence_to_edit, sequence_right = strip_saturation_parentheses(target_sequence)
    sm_target_sequence_list = []
    sm_target_name_list = []
    sm_output_target_name_list = []
    sm_allowed_range_list = []
    allowed_range = (len(sequence_left), len(sequence_left) + len(sequence_to_edit))

    def append_sm_target(output_target_name, sm_target_sequence, internal_suffix=None):
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
        sm_allowed_range_list.append(allowed_range)

    if (sm_type == 'aa') or (sm_type == 'aa_subs'):
        for base_index in range(0, len(sequence_to_edit), 3):
            codon_ref = sequence_to_edit[base_index:base_index + 3]
            if len(codon_ref) == 3:
                aa_ref = codon_dict[codon_ref][1]
                inner_left = sequence_to_edit[:base_index]
                inner_right = sequence_to_edit[base_index + 3:]
                for aa_edit in [x for x in aa2codon if x != aa_ref]:
                    for codon_edit, _usage in sorted_codons_by_difference(codon_ref, aa2codon[aa_edit]):
                        output_name = '%s_%s_%sto%s' % (target_name, str(int(base_index / 3 + 1)), aa_ref, aa_edit)
                        append_sm_target(output_name, sequence_left + inner_left + '(%s/%s)' % (codon_ref, codon_edit) + inner_right + sequence_right, '%s_%s' % (aa_edit, codon_edit))

    if (sm_type == 'aa') or (sm_type == 'aa_ins'):
        for base_index in range(0, len(sequence_to_edit), 3):
            codon_ref = sequence_to_edit[base_index:base_index + 3]
            if len(codon_ref) == 3:
                aa_ref = codon_dict[codon_ref][1]
                inner_left = sequence_to_edit[:base_index]
                inner_right = sequence_to_edit[base_index + 3:]
                for aa_edit in [x for x in aa2codon if x != 'X']:
                    for codon_edit, _usage in sorted_codons_by_usage(aa2codon[aa_edit]):
                        output_name = '%s_%s_%sto%s%s' % (target_name, str(int(base_index / 3 + 1)), aa_ref, aa_ref, aa_edit)
                        append_sm_target(output_name, sequence_left + inner_left + codon_ref + '(+%s)' % codon_edit + inner_right + sequence_right, '%s_%s' % (aa_edit, codon_edit))

    if (sm_type == 'aa') or (sm_type == 'aa_dels'):
        for base_index in range(0, len(sequence_to_edit), 3):
            codon_ref = sequence_to_edit[base_index:base_index + 3]
            seq_to_edit_and_right = sequence_to_edit + sequence_right
            codon2_ref = seq_to_edit_and_right[base_index + 3:base_index + 6]
            if len(codon_ref) == 3 and len(codon2_ref) == 3:
                aa_ref = codon_dict[codon_ref][1]
                aa2_ref = codon_dict[codon2_ref][1]
                inner_left = sequence_to_edit[:base_index]
                inner_right = sequence_to_edit[base_index + 3:]
                output_name = '%s_%s_%s%sto%s' % (target_name, str(int(base_index / 3 + 1)), aa_ref, aa2_ref, aa2_ref)
                append_sm_target(output_name, sequence_left + inner_left + '(-%s)' % codon_ref + inner_right + sequence_right, codon_ref)

    if sm_type == 'aa_silent':
        for base_index in range(0, len(sequence_to_edit), 3):
            codon_ref = sequence_to_edit[base_index:base_index + 3]
            if len(codon_ref) == 3:
                aa_ref = codon_dict[codon_ref][1]
                inner_left = sequence_to_edit[:base_index]
                inner_right = sequence_to_edit[base_index + 3:]
                codon_edit_list = [x for x in aa2codon[aa_ref] if x[0] != codon_ref]
                for codon_edit, _usage in sorted_codons_by_difference(codon_ref, codon_edit_list):
                    output_name = '%s_%s_%sto%s' % (target_name, str(int(base_index / 3 + 1)), aa_ref, aa_ref)
                    append_sm_target(output_name, sequence_left + inner_left + '(%s/%s)' % (codon_ref, codon_edit) + inner_right + sequence_right, codon_edit)

    if sm_type == 'base':
        for base_index in range(len(sequence_to_edit)):
            base_ref = sequence_to_edit[base_index]
            inner_left = sequence_to_edit[:base_index]
            inner_right = sequence_to_edit[base_index + 1:]
            for base_edit in [x for x in ['A', 'T', 'C', 'G'] if x != base_ref.upper()]:
                output_name = '%s_%s_%sto%s' % (target_name, str(base_index + 1), base_ref, base_edit)
                append_sm_target(output_name, sequence_left + inner_left + '(%s/%s)' % (base_ref, base_edit) + inner_right + sequence_right)

    return sm_target_name_list, sm_target_sequence_list, sm_output_target_name_list, sm_allowed_range_list


def read_targets(path):
    targets = []
    with open(path, 'r', newline='') as handle:
        next(handle)
        if path.lower().endswith('.csv'):
            reader = csv.reader(handle)
            for row in reader:
                if not row:
                    continue
                if len(row) < 2:
                    logger.error('A row in %s does not have at least 2 columns.', path)
                    sys.exit(1)
                targets.append((row[0].strip(), row[1].strip().upper()))
        elif path.lower().endswith('.txt'):
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 2:
                    logger.error('Line %s in %s does not have at least 2 columns.', line, path)
                    sys.exit(1)
                targets.append((parts[0].strip(), parts[1].strip().upper()))
        else:
            logger.error('Input file %s does not end with .txt or .csv.', path)
            sys.exit(1)
    return targets


def build_designs():
    target_designs = []
    for target_name, target_sequence in read_targets(file_in):
        if saturation_mutagenesis:
            sm_names, sm_seqs, sm_output_names, sm_allowed_ranges = saturating_mutagenesis_input_sequences(target_name, target_sequence, saturation_mutagenesis)
            for sm_name, sm_seq, output_name, allowed_range in zip(sm_names, sm_seqs, sm_output_names, sm_allowed_ranges):
                processed = process_sequence(sm_seq)
                target_designs.append({
                    'target_name': output_name,
                    'internal_target_name': sm_name,
                    'target_sequence': sm_seq,
                    'saturation_mutagenesis': saturation_mutagenesis,
                    'silent_allowed_range': allowed_range,
                    **processed,
                })
        else:
            processed = process_sequence(target_sequence)
            target_designs.append({
                'target_name': target_name,
                'internal_target_name': target_name,
                'target_sequence': target_sequence,
                'saturation_mutagenesis': '',
                'silent_allowed_range': None,
                **processed,
            })
    return target_designs


def main():
    designs = build_designs()
    if len(designs) == 0:
        logger.error('Input file %s does not have any entries. Include a header: target_name,target_sequence.', file_in)
        sys.exit(1)

    out_file = os.path.join(out_dir, '%s_DMSDesign.csv' % str(time.strftime('%Y%m%d_%I.%M.%S', time.localtime())))
    fieldnames = [
        'Target_name', 'Internal_target_name', 'Target_sequence', 'Design_number', 'Edit_type',
        'Saturation_mutagenesis', 'Reference_sequence', 'Edit_sequence', 'Edit_sequence_with_silent_mutations',
        'Edit_start_in_ref', 'Edit_span_length_ref', 'Edit_span_length_edit',
        'Silent_mutations_requested', 'Silent_mutation_mode', 'Silent_mutations_added', 'Silent_mutation_annotation', 'Silent_mutation_positions',
        'Silent_mutation_allowed_region', 'Silent_barcode_id', 'Silent_barcode_pattern', 'Silent_barcode_reuse_count',
    ]
    with open(out_file, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for idx, design in enumerate(designs, start=1):
            if silent_mutation_mode == "barcode":
                (
                    silent_edit_sequence,
                    silent_annotation,
                    silent_positions,
                    silent_added,
                    barcode_pattern,
                    barcode_id,
                    barcode_reuse_count,
                ) = apply_silent_barcode(
                    design["edit_sequence"],
                    design["reference_sequence"],
                    silent_mutation,
                    barcode_id=idx - 1,
                    allowed_range=design["silent_allowed_range"],
                )
            else:
                silent_edit_sequence, silent_annotation, silent_positions, silent_added = add_silent_mutations(
                    design["edit_sequence"],
                    design["reference_sequence"],
                    silent_mutation,
                    allowed_range=design["silent_allowed_range"],
                    mode=silent_mutation_mode,
                )
                barcode_pattern = ""
                barcode_id = ""
                barcode_reuse_count = ""
            allowed_region = '' if design['silent_allowed_range'] is None else '%d-%d' % (design['silent_allowed_range'][0], design['silent_allowed_range'][1] - 1)
            writer.writerow({
                'Target_name': design['target_name'],
                'Internal_target_name': design['internal_target_name'],
                'Target_sequence': design['target_sequence'],
                'Design_number': idx,
                'Edit_type': edit_type_from_target(design['target_sequence']),
                'Saturation_mutagenesis': design['saturation_mutagenesis'],
                'Reference_sequence': design['reference_sequence'],
                'Edit_sequence': design['edit_sequence'],
                'Edit_sequence_with_silent_mutations': silent_edit_sequence,
                'Edit_start_in_ref': design['edit_start_in_ref'],
                'Edit_span_length_ref': design['edit_span_length_w_ref'],
                'Edit_span_length_edit': design['edit_span_length_w_edit'],
                'Silent_mutations_requested': silent_mutation,
                'Silent_mutation_mode': silent_mutation_mode,
                'Silent_mutations_added': silent_added,
                'Silent_mutation_annotation': silent_annotation,
                'Silent_mutation_positions': silent_positions,
                'Silent_mutation_allowed_region': allowed_region,
                'Silent_barcode_id': barcode_id,
                'Silent_barcode_pattern': barcode_pattern,
                'Silent_barcode_reuse_count': barcode_reuse_count,
            })
    logger.info('DMSDesign completed. Wrote %s designs to %s', len(designs), out_file)


if __name__ == '__main__':
    main()
