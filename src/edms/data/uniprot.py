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
import pandas as pd
import os
import urllib.request
import urllib.error
from ..utils import mkdir

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