"""
Utilities for working with PDB / mmCIF structures and computing residue–residue
neighbors using MDAnalysis.

This module provides:
- retrieve(): download PDB or mmCIF files from the RCSB PDB API.
- compute_residue_neighbors(): compute residue–residue neighbors (protein, DNA,
  RNA, ligands) using MDAnalysis, with optional confidence weighting (e.g.,
  for AlphaFold models using pLDDT stored in the B-factor field).

MDAnalysis supports multiple structure/topology formats including PDB and
mmCIF; mmCIF is generally preferred for modern PDB entries because it
Avoids many of the PDB format limitations (field truncation, large atom
counts, etc.).
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Literal, Optional, Union

import io
import os
import urllib.error
import urllib.request

import numpy as np
import pandas as pd

from ..utils import mkdir

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO

_PathLike = Union[str, os.PathLike]


def retrieve(id: str, suf: str=None, dir: str=None, config: bool=True, base_url: str="https://files.rcsb.org/download") -> Path:
    """
    retrieve(): Download a PDB or mmCIF file from the RCSB PDB REST API.

    Parameters:
    id (str): 4-character PDB accession (e.g. "1CRN", case-insensitive).
    suf (str, optional): Output file suffix (e.g., .pdb or .cif (PDBx/mmCIF)).
    dir (str, optional): Save file to specified directory.
    config (bool, optional): Save to configuration directory (Default: True).
    base_url (str, optional): Base URL for the RCSB PDB REST API download endpoint (Default: "https://files.rcsb.org/download").
    """
    # Determine file format from filename if provided
    if suf is None: suf = ".pdb"  # Default to PDB
    elif suf != ".cif" and suf != ".pdb": raise ValueError(f'suf must end with ".cif" or ".pdb", not {suf!r}')

    # Download PDBx/mmCIF or PDB from RCSB PDB REST API
    id = id.lower()
    url = f"{base_url.rstrip('/')}/{id}{suf}"
    try:
        with urllib.request.urlopen(url) as resp:
            # UniProt flat files are UTF-8 text
            text = resp.read().decode("utf-8")
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Failed to download PDB entry {id} "
                           f"(HTTP {e.code}) from {url}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"Failed to reach RCSB PDB at {url}: {e.reason}") from e
    
    # Save file to...
    file = f"{id}{suf}"
    if dir is not None: # specified directory
        mkdir(dir) # Ensure directory exists
        out_path = Path(dir) / file if dir else Path(file)
        out_path.write_text(text)

    if config==True: # config directory
        dir = os.path.expanduser("~/.config/edms/PDB")
        mkdir(dir) # Ensure directory exists
        out_path = Path(dir) / file
        out_path.write_text(text)


def _infer_confidence_per_residue(u: "mda.Universe", selection: str) -> np.ndarray:
    """
    Infer a per-residue confidence value from B-factors for a given selection.

    This is useful for AlphaFold models, where pLDDT is stored in the B-factor
    field for each atom. For each residue, we take the mean B-factor across its
    atoms.

    Parameters
    ----------
    u:
        MDAnalysis Universe.
    selection:
        Atom selection string (e.g. "protein").

    Returns
    -------
    np.ndarray
        One confidence value per residue in the selection (order matches
        u.select_atoms(selection).residues).

    Notes
    -----
    If B-factors are all zero or missing, this will return an array of ones.
    """
    ag = u.select_atoms(selection)
    residues = ag.residues
    conf = np.empty(len(residues), dtype=float)

    for i, res in enumerate(residues):
        bvals = res.atoms.tempfactors
        if bvals.size == 0:
            conf[i] = 1.0
        else:
            mean_b = float(np.mean(bvals))
            # If all zero (e.g. no confidence info), fall back to 1.0
            conf[i] = mean_b if mean_b > 0 else 1.0

    return conf


def compute_residue_neighbors(
    structure: _PathLike,
    selection: str = "protein",
    partner_selection: Optional[str] = None,
    cutoff: float = 0,
    weigh_by_confidence: bool = False,
    confidence_selection: Optional[str] = None,
) -> pd.DataFrame:
    """
    Compute residue–residue neighbors in a structure using MDAnalysis.

    This function computes pairwise distances between residues defined by
    `selection` (and, optionally, a distinct `partner_selection`) and returns
    all residue–residue pairs within a given distance cutoff. Residue positions
    are defined as the center of geometry of all atoms in each residue.

    Parameters
    ----------
    structure:
        Path to a structure file (PDB, PDBx/mmCIF, etc.) readable by MDAnalysis.
    selection:
        MDAnalysis atom selection string for the "query" residues
        (default: "protein").
    partner_selection:
        Optional atom selection string for the "partner" residues. If None,
        neighbors are computed within the same selection (all-vs-all protein
        neighbors, for example). To study protein–DNA interactions, you might
        use `selection="protein", partner_selection="nucleic"`.
    cutoff:
        Distance cutoff in Å for defining neighbors (default: 6 Å for protein–protein and 10 Å for protein–nucleic acid).
    weigh_by_confidence:
        If True, adjust distances by per-residue confidence values inferred
        from B-factors (e.g., AlphaFold pLDDT scores). Effective distance is
        defined as:

            d_eff = d / (conf_i * conf_j)

        where conf_i and conf_j are per-residue confidence values normalized
        to [0, 1] by dividing by 100 (for AlphaFold pLDDT-style values).
    confidence_selection:
        Optional selection string to use when estimating confidence values.
        If None, this defaults to `selection`. Only relevant when
        `weigh_by_confidence=True`.

    Returns
    -------
    pandas.DataFrame
        A tidy DataFrame with one row per residue–residue neighbor pair, with
        the following columns:

        - query_index: integer index of the query residue (0-based within selection)
        - query_resid: residue id (as in the structure)
        - query_resname: residue name (3-letter code)
        - partner_index: integer index of the partner residue
        - partner_resid: partner residue id
        - partner_resname: partner residue name
        - distance: center-of-geometry distance in Å
        - effective_distance: distance weighted by confidence (if enabled;
          otherwise equal to `distance`)

    Notes
    -----
    - Requires MDAnalysis to be installed. If MDAnalysis is not available,
      this function will raise an ImportError.
    - Works for protein–protein, protein–DNA/RNA, or any other selection that
      MDAnalysis understands.
    """
    if structure.endswith(".cif") or structure.endswith(".mmcif"): # Convert PDBx/mmCIF to PDB if needed
        message = f"Converting PDBx/mmCIF to PDB format for MDAnalysis...\nPDBx/mmCIF: {structure}"
        
        PDBIO().set_structure(MMCIFParser().get_structure("id", structure)) # Load PDBx/mmCIF file and set structure for PDBIO
        structure = ".".join(structure.split(".")[:-1]) + ".pdb"
        PDBIO().save(structure) # Save as PDB file
        
        message += f"\nPDB file: {structure}"
        print(message)
    
    elif structure.endswith(".pdb") == False:
        raise ValueError(f"Unsupported file format: {structure}.\nPlease provide a PDB or PDBx/mmCIF file.")

    structure = Path(structure)

    # Load universe (topology and coordinates)
    u = mda.Universe(str(structure))

    # Query selection
    ag_query = u.select_atoms(selection)
    query_residues = ag_query.residues

    # Partner selection: default to same as query
    if partner_selection is None:
        ag_partner = ag_query
    else:
        ag_partner = u.select_atoms(partner_selection)

    partner_residues = ag_partner.residues

    # Compute center-of-geometry for each residue
    query_coords = np.array([res.atoms.center_of_geometry() for res in query_residues])
    partner_coords = np.array([res.atoms.center_of_geometry() for res in partner_residues])

    # Compute pairwise distances
    if partner_coords.shape[0] == 0 or query_coords.shape[0] == 0: # No atoms in one of the selections; return empty DataFrame
        print(f"Warning: No atoms found in one of the selections (query: {selection}, partner: {partner_selection}) for compute_residue_neighbors(); returning empty DataFrame.")
        return pd.DataFrame(
            columns=[
                "query_index",
                "query_resid",
                "query_resname",
                "query_chain",
                "partner_index",
                "partner_resid",
                "partner_resname",
                "partner_chain",
                "distance",
                "effective_distance",
            ]
        )
    dist = distance_array(query_coords, partner_coords)

    # Confidence weighting
    if weigh_by_confidence:
        if confidence_selection is None:
            confidence_selection = selection

        conf_vals = _infer_confidence_per_residue(u, confidence_selection)
        # Normalize to [0, 1] assuming values like pLDDT (0–100)
        conf_vals = np.clip(conf_vals / 100.0, 1e-6, 1.0)

        # Map confidence values to query/partner residues by index
        if len(conf_vals) != len(query_residues):
            # If selections differ, we conservatively set all confidence to 1.0
            conf_query = np.ones(len(query_residues), dtype=float)
            conf_partner = np.ones(len(partner_residues), dtype=float)
        else:
            conf_query = conf_vals
            conf_partner = conf_vals if partner_selection is None else np.ones(len(partner_residues), dtype=float)

        # Broadcast to full matrix
        conf_mat = np.outer(conf_query, conf_partner)
        eff_dist = dist / conf_mat
    else:
        eff_dist = dist

    # Collect neighbor pairs within cutoff
    if cutoff <= 0:
        if partner_selection is None or partner_selection == "protein":
            cutoff = 6.0  # Default for protein–protein
        elif partner_selection == "nucleic":
            cutoff = 10.0  # Default for protein–nucleic acid
    q_idx, p_idx = np.where(dist <= cutoff)
    rows = []
    for i, j in zip(q_idx, p_idx):
        q_res = query_residues[i]
        p_res = partner_residues[j]
        if q_res == p_res: continue  # Skip self-interactions
        rows.append(
            {
                "query_index": int(i),
                "query_resid": int(q_res.resid),
                "query_resname": str(q_res.resname),
                "query_chain": q_res.segid,
                "partner_index": int(j),
                "partner_resid": int(p_res.resid),
                "partner_resname": str(p_res.resname),
                "partner_chain": p_res.segid,
                "distance": float(dist[i, j]),
                "effective_distance": float(eff_dist[i, j]),
            }
        )

    return pd.DataFrame(rows)

BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

def is_backbone_atom(atom) -> bool:
    return atom.name in BACKBONE_ATOMS

def compute_residue_contacts(
    structure: _PathLike,
    selection: str = "protein",
    partner_selection: Optional[str] = None,
    cutoff: float = 3.5,
    include_hydrogens: bool = False,
    exclude_backbone_backbone: bool = True,
) -> pd.DataFrame:
    """
    Compute residue–residue contacts using minimum atom–atom heavy-atom distances.

    Parameters
    ----------
    structure : path-like
        PDB file (MDAnalysis-readable). CIF is not supported unless converted to PDB.
    selection : str
        Atom selection for query residues (e.g., "protein").
    partner_selection : str or None
        Selection for partner residues. If None, uses the same as `selection`.
    cutoff : float
        Maximum heavy-atom distance (Å) for calling a contact.
    include_hydrogens : bool
        If True, include hydrogens. Default: heavy atoms only.
    exclude_backbone_backbone : bool
        If True, exclude backbone–backbone atom pairs from contact consideration.

    Returns
    -------
    DataFrame with columns:
        query_index, query_resid, query_resname, query_chain,
        partner_index, partner_resid, partner_resname, partner_chain,
        min_atom_distance
    """
    if structure.endswith(".cif") or structure.endswith(".mmcif"): # Convert PDBx/mmCIF to PDB if needed
        message = f"Converting PDBx/mmCIF to PDB format for MDAnalysis...\nPDBx/mmCIF: {structure}"
        
        PDBIO().set_structure(MMCIFParser().get_structure("id", structure)) # Load PDBx/mmCIF file and set structure for PDBIO
        structure = ".".join(structure.split(".")[:-1]) + ".pdb"
        PDBIO().save(structure) # Save as PDB file
        
        message += f"\nPDB file: {structure}"
        print(message)
    
    elif structure.endswith(".pdb") == False:
        raise ValueError(f"Unsupported file format: {structure}.\nPlease provide a PDB or PDBx/mmCIF file.")
    
    u = mda.Universe(str(structure))

    # Atom selections
    atom_sel = "all" if include_hydrogens else "not name H*"

    ag_query = u.select_atoms(f"{selection} and {atom_sel}")
    ag_partner = u.select_atoms(
        f"{partner_selection} and {atom_sel}"
        if partner_selection else f"{selection} and {atom_sel}"
    )

    query_residues = ag_query.residues
    partner_residues = ag_partner.residues

    rows = []

    for i, q_res in enumerate(query_residues):
        q_atoms = q_res.atoms
        q_coords = q_atoms.positions

        for j, p_res in enumerate(partner_residues):

            # Skip self-contact only if selections are identical
            if (partner_selection is None or partner_selection == selection) and i == j:
                continue
            # Skip sequential peptide neighbors (same chain)
            elif q_res.segid == p_res.segid and abs(q_res.resid - p_res.resid) == 1:
                continue

            p_atoms = p_res.atoms
            p_coords = p_atoms.positions

            # Determine minimum atom–atom distance
            dist_mat = distance_array(q_coords, p_coords)
            
            # Build mask of allowed (q_atom, p_atom) pairs
            allowed = np.ones(dist_mat.shape, dtype=bool)

            # Exclude backbone–backbone atom pairs if requested
            if exclude_backbone_backbone:
                # vectorized backbone mask (no Python loops)
                q_bb_mask = np.array([is_backbone_atom(a) for a in q_atoms])   # (nq,)
                p_bb_mask = np.array([is_backbone_atom(a) for a in p_atoms])   # (np,)
                bb_mask = np.outer(q_bb_mask, p_bb_mask)                       # (nq, np)

                # Exclude backbone–backbone pairs
                allowed &= ~bb_mask

            # If nothing is allowed, skip this residue pair
            if not allowed.any():
                continue

            # Mask disallowed pairs by setting their distance to +inf
            masked_dist = np.where(allowed, dist_mat, np.inf)

            # Find the minimum over allowed pairs
            min_flat = masked_dist.argmin()
            min_dist = masked_dist.flat[min_flat]

            # If even the best allowed pair is above cutoff, no contact
            if not np.isfinite(min_dist) or min_dist > cutoff:
                continue
            
            # Decode indices back to atom indices
            q_idx, p_idx = divmod(min_flat, masked_dist.shape[1])

            q_atom = q_atoms[q_idx]
            p_atom = p_atoms[p_idx]
                
            rows.append(
                {
                    "query_index": i,
                    "query_resid": int(q_res.resid),
                    "query_resname": q_res.resname,
                    "query_chain": q_res.segid,
                    "query_atom": q_atom.name,

                    "partner_index": j,
                    "partner_resid": int(p_res.resid),
                    "partner_resname": p_res.resname,
                    "partner_chain": p_res.segid,
                    "partner_atom": p_atom.name,

                    "min_atom_distance": float(min_dist),
                }
            )

    return pd.DataFrame(rows)