'''
Module: plate.py
Author: Marc Zepeda
Created: 2026-01-04 
Description: Plates for biological experiments.

Usage:
[Helpers]
- WELL_SHAPES: Dict[int, Tuple[int, int]]
- _is_blank_row(row: List[str]) -> bool
- _looks_like_plate_header(row: List[str]) -> bool
- _normalize_cell(x: str) -> str
- _infer_plate_size(n_rows_present: int, n_cols_present: int) -> Optional
- PlateBlock (dataclass)

[CSV Parser]
- parse_csv_to_tidy(): parse a CSV that contains one or more plate blocks into tidy format

'''
# Imports
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Dict, Any
from ..utils import mkdir

import pandas as pd

# Helpers
WELL_SHAPES: Dict[int, Tuple[int, int]] = {
    6: (2, 3),
    12: (3, 4),
    24: (4, 6),
    48: (6, 8),
    96: (8, 12),
}

def _is_blank_row(row: List[str]) -> bool:
    return all((c or "").strip() == "" for c in row)

def _looks_like_plate_header(row: List[str]) -> bool:
    """
    Heuristic: first cell non-empty and at least one numeric-ish column header in the rest of the row.
    Examples: ["96-well 1","1","2",...]
    """
    if not row or (row[0] or "").strip() == "":
        return False
    # any header cell (after first) looks like an integer
    for c in row[1:]:
        if re.fullmatch(r"\s*\d+\s*", (c or "")):
            return True
    # also accept if it literally says "well" in top-left and has any non-empty later
    if re.search(r"well", row[0], re.IGNORECASE) and any((c or "").strip() for c in row[1:]):
        return True
    return False

def _normalize_cell(x: str) -> str:
    return (x or "").strip()

def _infer_plate_size(n_rows_present: int, n_cols_present: int) -> Optional[int]:
    # returns 6/12/24/48/96 if it matches known shape
    for wells, (r, c) in WELL_SHAPES.items():
        if (n_rows_present, n_cols_present) == (r, c):
            return wells
    return None

@dataclass
class PlateBlock:
    name: str
    header_cols: List[str]          # column labels aligned to positions 1..N
    rows: List[Tuple[str, List[str]]]  # (row_label, values aligned to positions 1..N)

def parse_csv_to_tidy(
    csv_path: str | Path,
    default_plate_name: str = "{wells}-well {index}",
    skip_empty_cells: bool = True,
    replicate_requires_value: bool = True,
    out_path: str = None,
) -> pd.DataFrame:
    """
    parse_csv_to_tidy(): parse a CSV that contains one or more plate blocks into tidy format.
    """
    csv_path = Path(csv_path)

    # Read raw CSV into a rectangular-ish list of lists
    with csv_path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        grid = [list(r) for r in reader]

    # Normalize row lengths
    max_len = max((len(r) for r in grid), default=0)
    for r in grid:
        if len(r) < max_len:
            r.extend([""] * (max_len - len(r)))

    # Find plate header row indices
    header_idxs: List[int] = []
    for i, row in enumerate(grid):
        if _looks_like_plate_header(row):
            header_idxs.append(i)

    if not header_idxs:
        raise ValueError("No plate header rows detected. Expect a row like: ['96-well 1','1','2',...]")

    # Add sentinel end
    header_idxs_sorted = sorted(header_idxs)
    header_idxs_sorted.append(len(grid))

    blocks: List[PlateBlock] = []

    for b in range(len(header_idxs_sorted) - 1):
        h = header_idxs_sorted[b]
        end = header_idxs_sorted[b + 1]

        header_row = [_normalize_cell(x) for x in grid[h]]
        raw_name = header_row[0]

        # Collect data rows until blank row OR next header
        data_rows: List[List[str]] = []
        for r in range(h + 1, end):
            if _is_blank_row(grid[r]):
                break
            if _looks_like_plate_header(grid[r]):
                break
            data_rows.append([_normalize_cell(x) for x in grid[r]])

        if not data_rows:
            continue

        # Determine max column extent from actual data (don’t invent missing cols)
        max_data_col = 0
        for dr in data_rows:
            for j in range(1, len(dr)):
                if dr[j] != "":
                    max_data_col = max(max_data_col, j)
        if max_data_col == 0:
            continue

        # Build header col labels up to max_data_col
        header_cols: List[str] = []
        for j in range(1, max_data_col + 1):
            lbl = header_row[j] if j < len(header_row) and header_row[j] != "" else str(j)
            header_cols.append(lbl)

        # Trim each data row to max_data_col
        rows_struct: List[Tuple[str, List[str]]] = []
        for dr in data_rows:
            if not dr or dr[0] == "":
                continue
            row_lbl = dr[0]
            vals = dr[1 : max_data_col + 1]
            if len(vals) < max_data_col:
                vals = vals + ([""] * (max_data_col - len(vals)))
            rows_struct.append((row_lbl, vals))

        wells = _infer_plate_size(len(rows_struct), len(header_cols))
        if raw_name == "":
            raw_name = default_plate_name.format(wells=wells or "unknown", index=(len(blocks) + 1))

        blocks.append(PlateBlock(name=raw_name, header_cols=header_cols, rows=rows_struct))

    # Build tidy records (NO plate-based replicate logic here)
    records: List[Dict[str, Any]] = []
    for block in blocks:
        for row_lbl, vals in block.rows:
            for j, val in enumerate(vals, start=1):
                if skip_empty_cells and val == "":
                    continue
                col_lbl = block.header_cols[j - 1]
                well = f"{row_lbl}{col_lbl}" if re.fullmatch(r"\d+", str(col_lbl)) else f"{row_lbl}{j}"
                records.append(
                    {
                        "plate": block.name,
                        "col": str(col_lbl),
                        "row": str(row_lbl),
                        "well": well,
                        "value": val,
                    }
                )

    df = pd.DataFrame.from_records(records)
    if df.empty:
        # keep expected columns
        return pd.DataFrame(columns=["plate", "col", "row", "well", "replicate", "value"])

    # Optional: stable ordering before replicate assignment
    df = df.sort_values(["plate", "row", "col"], kind="stable").reset_index(drop=True)

    # Replicate counting ONLY based on repeated values across whole file
    if replicate_requires_value:
        mask = df["value"].astype(str).str.strip().ne("")
        df["replicate"] = 1
        df.loc[mask, "replicate"] = df.loc[mask].groupby("value").cumcount() + 1
    else:
        df["replicate"] = df.groupby("value").cumcount() + 1

    if out_path is not None:
        mkdir('/'.join(out_path.split('/')[:-1]))
        df.to_csv(out_path, index=False)

    return df