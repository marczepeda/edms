''' 
Module: signature.py
Author: Marc Zepeda
Created: 2025-08-30

Usage:
[Dataclass]
- SNV: single nucleotide variant (position, reference, alternative)
- Indel: insertion/deletion (position, insertion, deletion length)
- Signature: SNV and Indel tuples

[Signature Literals]
- _ALLOWED_CTORS: Signature, SNV, and Indel constructors
- SignatureParseError(ValueError): Not a Signature object
- parse_signature_literal(): Convert string into a Signature object, if applicable

[Signature Functions]
- concat_gapped_from_aligned(): Construct continuous gapped strings for the aligned region only, using PairwiseAlignment.aligned coordinate blocks.
                                Works for global or local alignments. Leading/trailing unaligned sequence is omitted (local) or represented by gaps 
                                (if you add that explicitly).
- left_align_indels(): Left-align simple indels within repeat context on the reference
- signature_from_alignment(): Get Signature from reference-query sequence alignment
'''

# Import packages
from typing import Tuple, List, Dict, Any
from dataclasses import dataclass, is_dataclass
import ast

# Dataclass
@dataclass(frozen=True)
class SNV:
    pos: int   # 0-based position on WT
    ref: str
    alt: str

@dataclass(frozen=True)
class Indel:
    pos: int   # 0-based left-aligned ref position where the event occurs
    ins: str   # inserted sequence ("" if deletion)
    dellen: int  # deleted length (0 if insertion)

@dataclass(frozen=True)
class Signature:
    snvs: Tuple[SNV, ...]
    indels: Tuple[Indel, ...]

# Signature Literals
_ALLOWED_CTORS: Dict[str, type] = {
    "Signature": Signature,
    "SNV": SNV,
    "Indel": Indel,
}

class SignatureParseError(ValueError):
    pass

def parse_signature_literal(text: str) -> Signature:
    """
    parse_signature_literal(): Convert string into a Signature object, if applicable

    Parameters:
    text (str): text that might represent a Signature object
    """
    try:
        node = ast.parse(text, mode="eval").body
    except SyntaxError as e:
        return text   # not even valid Python expr → leave unchanged
        
    def conv(n: ast.AST) -> Any:
        # Signature(...), SNV(...), Indel(...)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name):
            name = n.func.id
            if name not in _ALLOWED_CTORS:
                raise SignatureParseError(f"Disallowed constructor: {name}")
            if n.args:
                # Keep it strict: only keyword args to avoid positional ambiguity
                raise SignatureParseError("Positional arguments are not allowed.")
            kwargs = {kw.arg: conv(kw.value) for kw in n.keywords}
            cls = _ALLOWED_CTORS[name]
            obj = cls(**kwargs)
            # Optional sanity check: ensure we created the right dataclass
            if not is_dataclass(obj):
                raise SignatureParseError(f"{name} is not a dataclass.")
            return obj

        # Tuples: ( ... )
        if isinstance(n, ast.Tuple):
            return tuple(conv(elt) for elt in n.elts)

        # Simple literals
        if isinstance(n, ast.Constant):
            if isinstance(n.value, (int, str)) or n.value is None:
                return n.value
            raise SignatureParseError(f"Disallowed constant: {n.value!r}")

        # Forbid everything else (Names, BinOps, f-strings, lists, dicts, etc.)
        raise SignatureParseError(f"Disallowed syntax: {type(n).__name__}")

    obj = conv(node)
    if not isinstance(obj, Signature):
        return text   # valid expr but not our form → leave unchanged
    return obj if isinstance(obj, Signature) else text

# Signature Functions
def concat_gapped_from_aligned(ref_seq: str, query_seq: str, alignment) -> Tuple[str, str]:
    """
    concat_gapped_from_aligned(): Construct continuous gapped strings for the aligned region only, using PairwiseAlignment.aligned coordinate blocks.
                                  Works for global or local alignments. Leading/trailing unaligned sequence is omitted (local) or represented by gaps 
                                  (if you add that explicitly).
    
    Parameters:
    ref_seq (str): reference sequence
    query_seq (str): query sequence
    alignment: PairwiseAligner alignment
    """
    ref_blocks = alignment.aligned[0]  # array of [start, end) on WT
    query_blocks = alignment.aligned[1] # array of [start, end) on read
    assert len(ref_blocks) == len(query_blocks), "Mismatched aligned block counts"

    ref_parts, query_parts = [], []

    for i, ((rs, re), (qs, qe)) in enumerate(zip(ref_blocks, query_blocks)):
        if i > 0:
            prs, pre = ref_blocks[i-1]
            pqs, pqe = query_blocks[i-1]

            # Gap(s) between previous and current blocks
            # Deletion relative to WT (gap in read)
            if rs > pre and qs == pqe:
                ref_parts.append(ref_seq[pre:rs])
                query_parts.append("-" * (rs - pre))

            # Insertion relative to WT (gap in reference)
            if qs > pqe and rs == pre:
                ref_parts.append("-" * (qs - pqe))
                query_parts.append(query_seq[pqe:qs])

            # If both advanced, that would indicate an unaligned diagonal jump
            # which shouldn't happen in a standard pairwise alignment.
            if rs > pre and qs > pqe:
                # Defensive: treat it as mismatches (diagonal) – append them aligned.
                span = min(rs - pre, qs - pqe)
                ref_parts.append(ref_seq[pre:pre+span])
                query_parts.append(query_seq[pqe:pqe+span])
                # Any remainder will be handled by the gap cases above
                if rs - pre > span:
                    ref_parts.append(ref_seq[pre+span:rs])
                    query_parts.append("-" * (rs - (pre + span)))
                if qs - pqe > span:
                    ref_parts.append("-" * (qs - (pqe + span)))
                    query_parts.append(query_seq[pqe+span:qs])

        # Now the aligned block itself (matches+mismatches)
        ref_parts.append(ref_seq[rs:re])
        query_parts.append(query_seq[qs:qe])

    ref_gapped = "".join(ref_parts)
    query_gapped = "".join(query_parts)
    return ref_gapped, query_gapped

def left_align_indels(indels: List[Indel], ref: str) -> List[Indel]:
    """
    left_align_indels(): Left-align simple indels within repeat context on the reference
    
    Parameters:
    indels (list[Indel]): list of Indel objects
    ref (str): reference sequence
    """
    out = []
    for ind in indels:
        if ind.dellen > 0:  # deletion
            start = ind.pos
            # shift left while previous base equals the rightmost deleted base
            del_seq = ref[ind.pos:ind.pos+ind.dellen]
            while start > 0 and ind.pos > 0 and ref[start-1] == del_seq[-1]:
                start -= 1
            out.append(Indel(pos=start, ins="", dellen=ind.dellen))
        elif ind.ins:  # insertion
            start = ind.pos
            # shift left while previous base equals the last base of insertion
            while start > 0 and ref[start-1] == ind.ins[-1]:
                start -= 1
            out.append(Indel(pos=start, ins=ind.ins, dellen=0))
        else:
            out.append(ind)
    
    # merge adjacent identical events if any
    out.sort(key=lambda x: (x.pos, x.dellen, x.ins))
    merged = []
    for ind in out:
        if merged and ind == merged[-1]:
            continue
        merged.append(ind)
    
    return merged

def signature_from_alignment(alignment, ref_seq: str, query_seq: str) -> Signature:
    '''
    signature_from_alignment(): Get Signature from reference-query sequence alignment

    Parameters:
    alignment: PairwiseAligner Alignment
    ref_seq (str): reference sequence
    query_seq (str): query sequence
    '''
    ref_g, query_g = concat_gapped_from_aligned(alignment=alignment,query_seq=query_seq,ref_seq=ref_seq)

    snvs: List[SNV] = []
    indels: List[Indel] = []

    ref_pos = -1  # will increment before use, so start at -1
    for r_base, q_base in zip(ref_g, query_g):
        if r_base != "-":
            ref_pos += 1

        if r_base == "-" and q_base != "-":
            # insertion relative to WT at current ref_pos (between ref_pos and ref_pos+1)
            # accumulate contiguous insertion
            ins_seq = q_base
            # continue through subsequent gap run
            # (we're iterating pairwise, so we can't lookahead easily; instead, we’ll compress later)
            indels.append(Indel(pos=ref_pos+1, ins=ins_seq, dellen=0))
        elif r_base != "-" and q_base == "-":
            # deletion relative to WT starting at this ref_pos
            # count contiguous deletion length
            indels.append(Indel(pos=ref_pos, ins="", dellen=1))
        elif r_base != "-" and q_base != "-" and r_base.upper() != q_base.upper():
            snvs.append(SNV(pos=ref_pos, ref=r_base.upper(), alt=q_base.upper()))
        # matches ignored

    # compress adjacent identical indel steps into single events
    compressed: Dict[Tuple[int, str, int], Indel] = {}
    for ind in indels:
        key = (ind.pos, ind.ins, int(ind.dellen > 0))
        if key in compressed:
            prev = compressed[key]
            if ind.dellen > 0:
                # extend deletion length if same start
                if ind.pos == prev.pos:
                    compressed[key] = Indel(pos=prev.pos, ins="", dellen=prev.dellen + 1)
                else:
                    compressed[(ind.pos, ind.ins, 1)] = ind
            else:
                # insertion pieces at same position -> concatenate
                compressed[key] = Indel(pos=prev.pos, ins=prev.ins + ind.ins, dellen=0)
        else:
            compressed[key] = ind
    indels_list = list(compressed.values())

    # left-align in repeat contexts for stability
    indels_list = left_align_indels(indels_list, ref_seq)

    snvs.sort(key=lambda s: s.pos)
    indels_list.sort(key=lambda d: (d.pos, d.dellen, d.ins))
    return Signature(snvs=tuple(snvs), indels=tuple(indels_list))