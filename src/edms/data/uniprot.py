"""
Module: uniprot.py
Author: Marc Zepeda
Created: 2024-12-01
Description: Utilities for parsing UniProt Feature viewer JSON exports.

This module is designed around (1) the JSON structure returned by the
UniProt "Feature viewer" export (e.g. P55317.json) and (2) the UniProtKB flat file text format (e.g. P55317.txt) in order to retrieve feature (FT) information.

Typical top-level keys (based on observed schema):
- entryType: str
- primaryAccession: str
- features: list[dict]
- extraAttributes: dict

Each feature usually has keys such as:
- type: str
- location: {
    "start": {"value": int, "modifier": str},
    "end":   {"value": int, "modifier": str},
  }
- description: str | None
- featureId: str | None
- evidences: list[{
    "evidenceCode": str,
    "source": str,
    "id": str,
  }]
- featureCrossReferences: list[{
    "database": str,
    "id": str,
  }]
- alternativeSequence: {
    "originalSequence": str | None,
    "alternativeSequences": list[str],
  }

The helpers below make it easy to:
- Load a UniProt Feature viewer JSON file.
- Work with a typed UniProtEntry/Feature representation in Python.
- Convert features to a tidy pandas.DataFrame for downstream analysis.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle, Circle, Arc
from matplotlib.patches import FancyArrowPatch
import os
import urllib.request
import urllib.error
from ..utils import mkdir
from ..gen import plot as p

# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

@dataclass
class Location:
    """Simple representation of a UniProt feature location.

    Attributes
    ----------
    start, end:
        1-based inclusive coordinates in the UniProt sequence.
    start_modifier, end_modifier:
        UniProt modifiers (e.g., "EXACT", "LESS_THAN"), if present.
    """

    start: int
    end: int
    start_modifier: Optional[str] = None
    end_modifier: Optional[str] = None


@dataclass
class Evidence:
    """Evidence supporting a UniProt feature annotation."""

    evidence_code: Optional[str] = None
    source: Optional[str] = None
    id: Optional[str] = None


@dataclass
class FeatureCrossReference:
    """Cross-references attached to a feature (e.g. dbSNP, PDB)."""

    database: str
    id: str


@dataclass
class AlternativeSequence:
    """Representation of the `alternativeSequence` subobject for variants."""

    original: Optional[str] = None
    alternatives: List[str] = field(default_factory=list)


@dataclass
class Feature:
    """A single UniProt feature entry."""

    type: str
    location: Location
    description: Optional[str] = None
    feature_id: Optional[str] = None
    evidences: List[Evidence] = field(default_factory=list)
    cross_references: List[FeatureCrossReference] = field(default_factory=list)
    alternative_sequence: Optional[AlternativeSequence] = None
    raw: Dict[str, Any] = field(default_factory=dict)

    @property
    def start(self) -> int:
        return self.location.start

    @property
    def end(self) -> int:
        return self.location.end

    @property
    def length(self) -> int:
        return self.location.end - self.location.start + 1


@dataclass
class UniProtEntry:
    """Top-level UniProt feature viewer entry."""

    accession: str
    entry_type: Optional[str] = None
    features: List[Feature] = field(default_factory=list)
    extra_attributes: Dict[str, Any] = field(default_factory=dict)
    raw: Dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Convenience methods
    # ------------------------------------------------------------------

    def filter_features(self, feature_type: str) -> List[Feature]:
        """Return features whose `type` matches `feature_type`.

        Parameters
        ----------
        feature_type:
            Case-sensitive feature type string as found in the JSON
            (e.g. "Natural variant", "Modified residue", "Region").
        """

        return [f for f in self.features if f.type == feature_type]

    def to_dataframe(self) -> "pd.DataFrame":  # type: ignore
        """Convert all features to a tidy pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
            One row per feature, with columns that are convenient for
            downstream analysis.

        Notes
        -----
        Requires pandas. If pandas is not installed, an ImportError is
        raised.
        """

        records: List[Dict[str, Any]] = []
        for feat in self.features:
            alt = feat.alternative_sequence
            evid_codes = [ev.evidence_code for ev in feat.evidences if ev.evidence_code]
            evid_sources = [ev.source for ev in feat.evidences if ev.source]
            evid_ids = [ev.id for ev in feat.evidences if ev.id]
            xref_db = [xr.database for xr in feat.cross_references]
            xref_id = [xr.id for xr in feat.cross_references]

            records.append(
                {
                    "accession": self.accession,
                    "entry_type": self.entry_type,
                    "feature_type": feat.type,
                    "feature_id": feat.feature_id,
                    "start": feat.start,
                    "end": feat.end,
                    "length": feat.length,
                    "description": feat.description,
                    "has_alternative_sequence": alt is not None,
                    "original_sequence": alt.original if alt else None,
                    "alternative_sequences": ",".join(alt.alternatives) if alt else None,
                    "evidence_codes": ";".join(evid_codes) if evid_codes else None,
                    "evidence_sources": ";".join(evid_sources) if evid_sources else None,
                    "evidence_ids": ";".join(evid_ids) if evid_ids else None,
                    "xref_databases": ";".join(xref_db) if xref_db else None,
                    "xref_ids": ";".join(xref_id) if xref_id else None,
                }
            )

        return pd.DataFrame.from_records(records)

# ---------------------------------------------------------------------------
# Retrieve flat file
# ---------------------------------------------------------------------------
def retrieve(accession: str, dir: str=None, config: bool=True, base_url: str = "https://rest.uniprot.org/uniprotkb") -> None:
    """
    retrieve(): Download a UniProt flat file entry via the REST API and save it to config and/or specified directory.

    Parameters:
    accession (str): UniProt accession (e.g. "P55317").
    dir (str, optional): Save file to specified directory.
    config (bool, optional): Save to configuration directory (Default: True).
    base_url (str, optional): Base URL for the UniProtKB REST API flat-file endpoint (Default: "https://rest.uniprot.org/uniprotkb").
    """
    # Download from UniProtKB REST API flat-file endpoint
    url = f"{base_url.rstrip('/')}/{accession}.txt"
    try:
        with urllib.request.urlopen(url) as resp:
            # UniProt flat files are UTF-8 text
            text = resp.read().decode("utf-8")
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Failed to download UniProt flat file for {accession} "
                           f"(HTTP {e.code}) from {url}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"Failed to reach UniProt REST API at {url}: {e.reason}") from e

    file = f"{accession}.txt" # Create filename

    # Save file to...
    if dir is not None: # specified directory
        mkdir(dir) # Ensure directory exists
        out_path = Path(dir) / file if dir else Path(file)
        out_path.write_text(text)

    if config==True: # config directory
        dir = os.path.expanduser("~/.config/edms/UniProt")
        mkdir(dir) # Ensure directory exists
        out_path = Path(dir) / file
        out_path.write_text(text)

# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

_JSONLike = Union[str, bytes, Dict[str, Any]]
_PathLike = Union[str, Path]

def _parse_flat_text(text: str) -> UniProtEntry:
    """Parse a UniProt flat file entry (text) into a UniProtEntry.

    This expects the standard UniProtKB/Swiss-Prot flat file format with
    ID/AC/FT/SQ lines, such as the P55317.txt example.
    """

    lines = text.splitlines()

    accession: Optional[str] = None
    entry_type: Optional[str] = None

    # ------------------------------------------------------------------
    # Parse accession and entry type from ID/AC lines
    # ------------------------------------------------------------------
    for line in lines:
        if line.startswith("ID"):
            # Example: "ID   FOXA1_HUMAN             Reviewed;         472 AA."
            if "Reviewed;" in line:
                entry_type = "Reviewed"
            elif "Unreviewed;" in line:
                entry_type = "Unreviewed"
        elif line.startswith("AC") and accession is None:
            # Example: "AC   P55317; B2R9H6; B7ZAP5; Q9H2A0;"
            rest = line[5:].strip()
            if rest:
                primary = rest.split(";")[0].strip()
                if primary:
                    accession = primary
        if accession is not None and entry_type is not None:
            break

    if accession is None:
        raise ValueError("Flat file does not contain an 'AC' line with an accession.")

    # ------------------------------------------------------------------
    # Parse FT feature table
    # ------------------------------------------------------------------
    features_raw: List[Dict[str, Any]] = []
    current_feature: Optional[Dict[str, Any]] = None
    current_qual_key: Optional[str] = None
    current_qual_value_parts: List[str] = []

    in_multiline_qual = False

    def flush_multiline_qual() -> None:
        nonlocal in_multiline_qual, current_qual_key, current_qual_value_parts
        if current_feature is not None and current_qual_key is not None:
            value = " ".join(part.strip() for part in current_qual_value_parts).strip()
            # Strip surrounding quotes if present
            if value.startswith('"') and value.endswith('"') and len(value) >= 2:
                value = value[1:-1]
            qualifiers = current_feature.setdefault("qualifiers", {})
            qualifiers[current_qual_key] = value
        in_multiline_qual = False
        current_qual_key = None
        current_qual_value_parts = []

    for line in lines:
        if not line.startswith("FT"):
            # When leaving FT section or hitting non-FT lines, flush any
            # unfinished multiline qualifier.
            if in_multiline_qual:
                flush_multiline_qual()
            continue

        # Slice after "FT   " (2 chars + 3 spaces). This is compatible
        # with standard UniProt formatting.
        content = line[5:]
        if not content.strip():
            continue

        # If this is the start of a new feature: non-space in the first column
        if content[0] != " ":
            # Finish previous qualifier if needed
            if in_multiline_qual:
                flush_multiline_qual()

            # Flush previous feature
            if current_feature is not None:
                features_raw.append(current_feature)

            # New feature line, columns ~6-13: key, then location
            key = content[:8].strip()
            loc_str = content[8:].strip()

            # Parse location such as "1..472" or "307"
            start: Optional[int] = None
            end: Optional[int] = None
            if ".." in loc_str:
                start_str, end_str = loc_str.split("..", 1)
                start_str = start_str.strip()
                end_str = end_str.strip()
                if start_str.isdigit():
                    start = int(start_str)
                # end may contain non-numeric suffix; extract leading digits
                end_digits = "".join(ch for ch in end_str if ch.isdigit())
                if end_digits:
                    end = int(end_digits)
            else:
                digits = "".join(ch for ch in loc_str if ch.isdigit())
                if digits:
                    start = end = int(digits)

            if start is None or end is None:
                # Skip malformed feature
                current_feature = None
                continue

            current_feature = {
                "type": key,
                "location": {"start": {"value": start}, "end": {"value": end}},
                "qualifiers": {},
            }
            current_qual_key = None
            current_qual_value_parts = []
            in_multiline_qual = False
            continue

        # From here on, we are in a continuation / qualifier line
        stripped = content.strip()

        # New qualifier
        if stripped.startswith("/") and "=" in stripped:
            # Flush previous multiline qualifier
            if in_multiline_qual:
                flush_multiline_qual()

            if current_feature is None:
                continue

            q = stripped[1:]
            qkey, qval = q.split("=", 1)
            qkey = qkey.strip()
            qval = qval.strip()

            # Start of quoted qualifier
            if qval.startswith('"') and not qval.endswith('"'):
                in_multiline_qual = True
                current_qual_key = qkey
                current_qual_value_parts = [qval]
            else:
                # Single-line qualifier
                if qval.startswith('"') and qval.endswith('"') and len(qval) >= 2:
                    qval = qval[1:-1]
                qualifiers = current_feature.setdefault("qualifiers", {})
                qualifiers[qkey] = qval
            continue

        # Continuation of a multiline qualifier
        if in_multiline_qual and current_qual_key is not None:
            current_qual_value_parts.append(stripped)
            # Check if this line closes the quote
            if stripped.endswith('"'):
                flush_multiline_qual()
            continue

        # Anything else in the FT section that doesn't match the above is
        # currently ignored.

    # Flush final multiline qualifier and feature
    if in_multiline_qual:
        flush_multiline_qual()
    if current_feature is not None:
        features_raw.append(current_feature)

    # Convert raw feature dicts into Feature dataclasses
    features: List[Feature] = []
    for f_raw in features_raw:
        loc = _parse_location(f_raw["location"])
        qualifiers: Dict[str, Any] = f_raw.get("qualifiers", {})
        note = qualifiers.get("note")
        feat_id = qualifiers.get("id")
        evidence_str = qualifiers.get("evidence")

        raw: Dict[str, Any] = {"qualifiers": qualifiers}
        if evidence_str is not None:
            raw["evidence_str"] = evidence_str

        features.append(
            Feature(
                type=f_raw.get("type", ""),
                location=loc,
                description=note,
                feature_id=feat_id,
                evidences=[],
                cross_references=[],
                alternative_sequence=None,
                raw=raw,
            )
        )

    return UniProtEntry(
        accession=accession,
        entry_type=entry_type,
        features=features,
        extra_attributes={},
        raw={"source": "flat", "accession": accession},
    )


def load_flat_file(obj: _PathLike) -> UniProtEntry:
    """Load a UniProt flat text entry (ID/AC/FT/SQ format) into a UniProtEntry.

    Parameters
    ----------
    obj:
        Path to a UniProt flat file (str or Path) or a raw string
        containing a single UniProt entry.
    """

    if isinstance(obj, (str, Path)):
        path = Path(str(obj))
        if path.exists():
            text = path.read_text()
        else:
            # Treat as raw text content
            text = str(obj)
    else:
        raise TypeError("obj must be a path or string containing a UniProt entry.")

    return _parse_flat_text(text)



def secondary_structure_from_flat_file(obj: _PathLike) -> pd.DataFrame:
    """
    Extract HELIX and STRAND features from a UniProt flat file into a DataFrame.

    Columns:
        accession: UniProt accession
        type: 'helix' or 'strand'
        start: 1-based start residue index
        end: 1-based end residue index
        name: 'helix #1', 'helix #2', ..., 'strand #1', ...

    Parameters
    ----------
    obj:
        Path to a UniProt flat file (e.g. 'P55317.txt') or a raw text
        containing a single UniProt entry.
    """
    # Reuse your existing loader
    entry = load_flat_file(obj)

    # Keep only HELIX and STRAND features
    ss_feats = [f for f in entry.features if f.type in {"HELIX", "STRAND"}]

    rows = []
    counters = {"HELIX": 0, "STRAND": 0}

    for feat in ss_feats:
        counters[feat.type] += 1
        idx = counters[feat.type]
        type_lower = "α-helix" if feat.type == "HELIX" else "β-strand"
        name = f"{type_lower} #{idx}"

        rows.append(
            {
                "accession": entry.accession,
                "type": type_lower,
                "start": feat.start,
                "end": feat.end,
                "name": name,
            }
        )

    return pd.DataFrame(rows)

def normalize_ptm_description(desc):
    if desc is None:
        return None
    desc = desc.lower()
    if "phospho" in desc:
        return "Phosphorylation"
    if "acetyl" in desc:
        return "Acetylation"
    if "methyl" in desc:
        return "Methylation"
    if "glycosylation" in desc:
        return "Glycosylation"
    if "ubiquitin" in desc:
        return "Ubiquitination"
    return desc

def ptms_from_flat_file(obj: _PathLike) -> pd.DataFrame:
    entry = load_flat_file(obj)

    ptm_feature_types = {"MOD_RES", "CARBOHYD", "LIPID", "DISULFID", "CROSSLNK", "INIT_MET"}

    rows = []
    for feat in entry.features:
        if feat.type in ptm_feature_types:
            rows.append({
                "accession": entry.accession,
                "type": feat.type,
                "start": feat.start,
                "end": feat.end,
                "description": feat.description,
                "normalized_description": normalize_ptm_description(feat.description),
            })

    return pd.DataFrame(rows)

def _parse_location(data: Dict[str, Any]) -> Location:
    """Parse a UniProt `location` object into a Location instance."""

    start_obj = data.get("start", {}) or {}
    end_obj = data.get("end", {}) or {}

    start_val = int(start_obj.get("value"))
    end_val = int(end_obj.get("value"))

    return Location(
        start=start_val,
        end=end_val,
        start_modifier=start_obj.get("modifier"),
        end_modifier=end_obj.get("modifier"),
    )


def _parse_evidences(data: Iterable[Dict[str, Any]] | None) -> List[Evidence]:
    """Parse a list of evidence objects."""

    if not data:
        return []

    out: List[Evidence] = []
    for ev in data:
        out.append(
            Evidence(
                evidence_code=ev.get("evidenceCode"),
                source=ev.get("source"),
                id=str(ev.get("id")) if ev.get("id") is not None else None,
            )
        )
    return out


def _parse_cross_references(
    data: Iterable[Dict[str, Any]] | None,
) -> List[FeatureCrossReference]:
    """Parse a list of featureCrossReferences objects."""

    if not data:
        return []

    out: List[FeatureCrossReference] = []
    for cr in data:
        db = cr.get("database")
        id_ = cr.get("id")
        if db is None or id_ is None:
            continue
        out.append(FeatureCrossReference(database=str(db), id=str(id_)))
    return out


def _parse_alt_sequence(data: Dict[str, Any] | None) -> Optional[AlternativeSequence]:
    """Parse the `alternativeSequence` subobject, if present."""

    if not data:
        return None

    original = data.get("originalSequence")
    alternatives_raw = data.get("alternativeSequences") or []
    alternatives = [str(a) for a in alternatives_raw]

    return AlternativeSequence(original=original, alternatives=alternatives)


def _parse_feature(data: Dict[str, Any]) -> Feature:
    """Parse a raw feature dictionary into a Feature dataclass."""

    loc = _parse_location(data.get("location", {}))

    evidences = _parse_evidences(data.get("evidences"))
    xrefs = _parse_cross_references(data.get("featureCrossReferences"))
    alt = _parse_alt_sequence(data.get("alternativeSequence"))

    return Feature(
        type=data.get("type", ""),
        location=loc,
        description=data.get("description"),
        feature_id=data.get("featureId"),
        evidences=evidences,
        cross_references=xrefs,
        alternative_sequence=alt,
        raw=data,
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_uniprot_feature_viewer_json(
    obj: _JSONLike,
) -> UniProtEntry:
    """Load a UniProt Feature viewer JSON export into a UniProtEntry.

    Parameters
    ----------
    obj:
        One of:
        - Path to a JSON file (str or Path).
        - Raw JSON string or bytes.
        - Already-parsed dict matching the UniProt Feature viewer schema.

    Returns
    -------
    UniProtEntry
        Parsed entry object with Features and metadata.
    """

    if isinstance(obj, (str, bytes)):
        # If this looks like a path that exists on disk, read from it;
        # otherwise treat it as raw JSON text.
        path = Path(str(obj))
        if path.exists():
            data = json.loads(path.read_text())
        else:
            data = json.loads(obj)
    elif isinstance(obj, dict):
        data = obj
    else:
        raise TypeError(
            "obj must be a path, JSON string, bytes, or dict; got " f"{type(obj)!r}"
        )

    entry_type = data.get("entryType")
    accession = data.get("primaryAccession") or data.get("accession")
    if accession is None:
        raise ValueError("JSON does not contain a 'primaryAccession' or 'accession' field.")

    features_raw: Sequence[Dict[str, Any]] = data.get("features", [])
    features = [_parse_feature(f) for f in features_raw]

    extra_attributes = data.get("extraAttributes") or {}

    return UniProtEntry(
        accession=str(accession),
        entry_type=entry_type,
        features=list(features),
        extra_attributes=extra_attributes,
        raw=data,
    )


def load_many_from_files(paths: Iterable[_PathLike]) -> List[UniProtEntry]:
    """Load multiple UniProt Feature viewer JSON files.

    Parameters
    ----------
    paths:
        Iterable of file paths.

    Returns
    -------
    list[UniProtEntry]
    """

    entries: List[UniProtEntry] = []
    for p in paths:
        entries.append(load_uniprot_feature_viewer_json(Path(p)))
    return entries


def entries_to_dataframe(entries: Iterable[UniProtEntry]) -> "pd.DataFrame":  # type: ignore
    """Convert multiple UniProtEntry objects to a single tidy DataFrame.

    Each feature from every entry becomes a row in the output.

    Notes
    -----
    Requires pandas. If pandas is not installed, an ImportError is raised.
    """

    frames = [e.to_dataframe() for e in entries]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


__all__ = [
    "Location",
    "Evidence",
    "FeatureCrossReference",
    "AlternativeSequence",
    "Feature",
    "UniProtEntry",
    "load_uniprot_feature_viewer_json",
    "load_many_from_files",
    "entries_to_dataframe",
    "load_flat_file",
]

# ---------------------------------------------------------------------------
# Secondary structure track
# ---------------------------------------------------------------------------

def draw_strand(
    ax: plt.Axes,
    start: float,
    end: float,
    y: float=0,
    height: float=25,
    color: str="gold",
    edgecolor: str="black",
    lw: float=2.0,
    head_frac: float=0.25,
    zorder: int=1,
):
    '''
    draw_strand(): Draw a strand on the axis.

    Parameters:
    ax : matplotlib axis
    start (float): start position in data coordinates
    end (float): end position in data coordinates
    y (float): y position in data coordinates
    height (float): height of the strand
    color (str): color of the strand
    edgecolor (str): color of the edge of the strand
    lw (float): line width of the strand
    head_frac (float): fraction of the strand length that is the head
    zorder (int): z-order of the strand
    '''
    length = end - start
    if length <= 0:
        return

    #head_length = 1
    head_length = max(length * head_frac, 5)
    tail_width = height
    head_width = height * 1.4

    arrow = FancyArrowPatch(
        (start, y),
        (end, y),
        arrowstyle=f"simple,tail_width={tail_width},"
                   f"head_width={head_width},"
                   f"head_length={head_length}",
        facecolor=color,
        edgecolor=edgecolor,
        linewidth=lw,
        mutation_scale=1,
        zorder=zorder,
    )
    ax.add_patch(arrow)

def draw_helix(
    ax: plt.Axes,
    start: float,
    end: float,
    y: float=0,
    height: float=0.6,
    color: str="turquoise",
    edgecolor: str="black",
    lw: float=2.0,
    left_cap_ls: str="-",
    right_front_ls: str="-",
    right_back_ls: str=":",
    right_back_alpha: float=0.8,
    zorder: int=1,
):
    '''
    draw_helix(): Draw a helix with a split right cap.
    
    Parameters:
    ax : matplotlib axis
    start (float): start position in data coordinates
    end (float): end position in data coordinates
    y (float): y position in data coordinates
    height (float): height of the helix
    color (str): color of the helix
    edgecolor (str): color of the edge of the helix
    lw (float): line width of the helix
    left_cap_ls (str): line style of the left cap
    right_front_ls (str): line style of the front right cap
    right_back_ls (str): line style of the back right cap
    right_back_alpha (float): alpha value of the back right cap
    zorder (int): z-order of the helix
    '''
    if end <= start:
        return

    length = end - start
    r = height / 2.0

    # very short helix: fallback blob
    if length <= height:
        circ = Circle(((start + end) / 2.0, y), radius=length / 2.0,
                      facecolor=color, edgecolor=edgecolor, linewidth=lw, zorder=zorder)
        ax.add_patch(circ)
        return

    # Fill capsule (no edge)
    rect = Rectangle((start + r, y - r), width=length - 2*r, height=height,
                     facecolor=color, edgecolor="none", linewidth=0, zorder=zorder)
    left = Circle((start + r, y), radius=r, facecolor=color, edgecolor="none", linewidth=0, zorder=zorder)
    right = Circle((end - r, y), radius=r, facecolor=color, edgecolor="none", linewidth=0, zorder=zorder)
    ax.add_patch(rect); ax.add_patch(left); ax.add_patch(right)

    # Body outline
    x0, x1 = start + r, end - r
    ax.plot([x0, x1], [y + r, y + r], color=edgecolor, lw=lw, ls='-', zorder=zorder + 0.1)
    ax.plot([x0, x1], [y - r, y - r], color=edgecolor, lw=lw, ls='-', zorder=zorder + 0.1)

    # Left cap outline (full)
    cxL = start + r
    ax.add_patch(Arc((cxL, y), width=2*r, height=2*r, angle=0,
                     theta1=-90, theta2=90, color=edgecolor, lw=lw, ls=left_cap_ls, zorder=zorder + 0.1))
    ax.add_patch(Arc((cxL, y), width=2*r, height=2*r, angle=0,
                     theta1=90, theta2=270, color=edgecolor, lw=lw, ls=left_cap_ls, zorder=zorder + 0.1))

    # Right cap outline split: right half solid, left half dotted
    cxR = end - r
    ax.add_patch(Arc((cxR, y), width=2*r, height=2*r, angle=0,
                     theta1=-90, theta2=90, color=edgecolor, lw=lw, ls=right_front_ls, zorder=zorder + 0.1))
    ax.add_patch(Arc((cxR, y), width=2*r, height=2*r, angle=0,
                     theta1=90, theta2=270, color=edgecolor, lw=lw, ls=right_back_ls,
                     alpha=right_back_alpha, zorder=zorder + 0.1))

def draw_loop(
    ax: plt.Axes,
    start: float,
    end: float,
    y: float=0,
    height: float=0.6,
    color: str="black",
    lw: float=2.0,
    cycles: int=None,
    resolution: int=1000,
    zorder: int=1,
):
    '''
    draw_loop(): Draw a loop as a sine wave.
    
    Parameters:
    ax : matplotlib axis
    start (float): start position in data coordinates
    end (float): end position in data coordinates
    y (float): y position in data coordinates
    height (float): height of the loop (peak-to-peak amplitude)
    color (str): color of the loop
    lw (float): line width of the loop
    cycles (int): number of sine wave cycles to draw; if None, automatically scaled based on length
    resolution (int): number of points to use for drawing the sine wave
    zorder (int): z-order of the loop
    '''
    if end <= start:
        return

    length = end - start
    if cycles is None:
        cycles = max(0.5, int(length / 2))

    x = np.linspace(start, end, resolution)
    amp = height / 2
    phase = 2 * np.pi * cycles * (x - start) / length
    y_wave = y + amp * np.sin(phase)

    ax.plot(x, y_wave, color=color, linewidth=lw,
            solid_capstyle="round", zorder=zorder)

def draw_ss_track(
    df: pd.DataFrame | str,
    y: float = 0,

    track_start: float | None = None,
    track_end: float | None = None,

    helix_height: float = 0.8,
    strand_height: float = 25,
    loop_height: float = 0.8,
    xpad: int = 5,

    helix_color: str = "turquoise",
    strand_color: str = "gold",
    loop_color: str = "black",
    edgecolor: str = "black",

    helix_lw: float = 2.0,
    strand_lw: float = 2.0,
    loop_lw: float = 2.0,

    loop_min_len: float = 1,

    base_zorder: int = 1,
    zorder_step: int = 1,

    helix_loop_pad: float = 0.0,
    strand_loop_pad: float = 0.25,

    figsize: tuple = (10, 2),
    dir: str | None = None,
    file: str | None = None,
    dpi: int = 0,
    show: bool = True
):
    """
    draw_ss_track(): Draw a secondary structure track

    Parameters:
    df (pd.DataFrame): UniProt accession (if saved to ~/.config/edms/UniProt) or file path for UniProt flat file. See edms.dat.uniprot.retrieve() or edms uniprot retrieve -h for more information.
    y (float): vertical position of the track
    track_start, track_end (float or None): optional start/end bounds for the track; if None, determined from data
    helix_height, strand_height, loop_height (float): heights for each element type; if any is None, defaults to `height`
    xpad (int): horizontal padding on the left and right of the track
    helix_color, strand_color, loop_color (str): colors for each element type
    edgecolor (str): color of the edges of helices and strands
    helix_lw, strand_lw, loop_lw (float): line widths for each element type
    loop_min_len (float): minimum length of loops to draw (for auto-drawn loops between elements)
    base_zorder (int): z-order for the first element; subsequent elements decrease by zorder_step to create a stacking effect
    helix_loop_pad, strand_loop_pad (float): padding to apply to the start/end of helices and strands when auto-drawing loops between elements, to make the sine wave connect more nicely; applied as a subtraction from the end of the current element and an addition to the start of the next element when determining loop start and end positions.
    figsize (tuple, optional): size of the figure to create
    dir (str, optional): output directory
    file (str, optional): output filename
    dpi (int, optional): figure dpi (Default: 1200)
    show (bool, optional): show plot (Default: True)
    """
    if isinstance(df, str):
        try: # from filepath
            df = secondary_structure_from_flat_file(obj=df)
        except:
            try: # from config
                for UniProt_file in os.listdir(os.path.expanduser('~/.config/edms/UniProt/')):
                    if df.lower() in UniProt_file.lower():
                        df = secondary_structure_from_flat_file(obj=f'{os.path.expanduser("~/.config/edms/UniProt")}/{UniProt_file}')
                        break
            except:
                raise FileNotFoundError(f"UniProt flat file not found: {df}.\nPlease provide a valid filename or UniProt accession (if saved to {os.path.expanduser('~/.config/edms/UniProt/')}) or file path for UniProt flat file. See edms.dat.uniprot.retrieve_flat_file() or edms uniprot retrieve -h for more information.")

    fig, ax = plt.subplots(figsize=figsize)

    df["start"] = df["start"].astype(float)
    df["end"] = df["end"].astype(float)
    df = df.sort_values("start").reset_index(drop=True)

    # Decide plotting bounds
    first_start = float(df.loc[0, "start"])
    last_end = float(df.loc[len(df) - 1, "end"])

    if track_start is None:
        track_start = first_start
    if track_end is None:
        track_end = last_end

    # ---- terminal loop: track_start -> first element ----
    if track_start < first_start:
        draw_loop(
            ax,
            start=track_start,
            end=first_start,
            y=y,
            height=loop_height,
            color=loop_color,
            lw=loop_lw,
            zorder=base_zorder + 1,   # keep terminal loop underneath everything
        )

    # ---- draw elements + internal loops ----
    for i, row in df.iterrows():
        current_z = base_zorder - i * zorder_step  # N->C decrease

        t = str(row["type"])
        s = float(row["start"])
        e = float(row["end"])

        is_helix = ("α" in t) or ("alpha" in t.lower()) or ("helix" in t.lower())
        is_strand = ("β" in t) or ("beta" in t.lower()) or ("strand" in t.lower()) or ("sheet" in t.lower())

        if is_helix:
            draw_helix(
                ax, s, e, y=y, height=helix_height,
                color=helix_color, edgecolor=edgecolor, lw=helix_lw, zorder=current_z
            )
        elif is_strand:
            draw_strand(
                ax, s, e, y=y, height=strand_height,
                color=strand_color, edgecolor=edgecolor, lw=strand_lw, zorder=current_z
            )
        else:
            draw_loop(
                ax, s, e, y=y, height=loop_height,
                color=loop_color, lw=loop_lw, zorder=current_z
            )

        # loop to next element
        if i < len(df) - 1:
            next_s = float(df.loc[i + 1, "start"])
            next_t = str(df.loc[i + 1, "type"])

            # Pad depending on type so the sine connects nicely
            pad_right = helix_loop_pad if is_helix else (strand_loop_pad if is_strand else 0.0)
            next_is_helix = ("α" in next_t) or ("alpha" in next_t.lower()) or ("helix" in next_t.lower())
            next_is_strand = ("β" in next_t) or ("beta" in next_t.lower()) or ("strand" in next_t.lower()) or ("sheet" in next_t.lower())
            pad_left = helix_loop_pad if next_is_helix else (strand_loop_pad if next_is_strand else 0.0)

            loop_s = e - pad_right
            loop_e = next_s + pad_left

            if loop_e > loop_s and (loop_e - loop_s) >= loop_min_len:
                # keep loop slightly underneath the "later" element so it doesn't cover it
                draw_loop(
                    ax, loop_s, loop_e, y=y, height=loop_height,
                    color=loop_color, lw=loop_lw, zorder=current_z
                )

    # ---- terminal loop: last element -> track_end ----
    if track_end > last_end:
        draw_loop(
            ax,
            start=last_end-.25,
            end=track_end,
            y=y,
            height=loop_height,
            color=loop_color,
            lw=loop_lw,
            zorder=base_zorder + len(df) * zorder_step - .5,
        )

    # Axis formatting
    xmin = min(track_start, first_start) - xpad
    xmax = max(track_end, last_end) + xpad

    # vertical limits: use helix height as the governing “track thickness”
    H = max(float(helix_height), float(loop_height))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-2,2)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_xlabel("Residue")

    p.save_fig(file=file, dir=dir, fig=ax.figure, dpi=dpi)
    if show:
        plt.show()