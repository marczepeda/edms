'''
Module: plates.py
Author: Marc Zepeda
Created: 2026-11-13
Description: Culture vessel formats and plate-aware transfection planning

A protocol written for a 96-well plate is the same protocol in a 6-well plate
with bigger numbers. This module holds the geometry of the vessels we actually
use, and wraps the transfection calculators in edms.bio.transfect so that
choosing a format rescales the DNA and the transfection volume instead of
leaving a stale 96-well number in the prose.

Scaling is by growth area, anchored on the format each calculator was written
for: edms.bio.transfect:PE3 on a 96-well plate, edms.bio.transfect:virus on a
6-well plate. The anchor amounts are unchanged, so a run in the anchor format
gets exactly the numbers it always did.

Not to be confused with edms.bio.plate, which parses plate maps into tidy
tables; this module is about the vessels themselves.

Usage:
[Formats]
- spec(): growth area, working media volume and seeding density for a format
- names(): canonical format names, smallest first
- table(): every format as a DataFrame
- scale(): growth-area ratio between two formats

[Planning]
- geometry(): one-row plan of the chosen format, for a protocol step
- PE3(): edms.bio.transfect:PE3 rescaled to a plate format
- virus(): edms.bio.transfect:virus rescaled to a plate format
'''
from __future__ import annotations

import pandas as pd

from ..bio import transfect

# Growth area (cm^2) and working media volume (mL) per well or dish.
# Areas are the usual Corning/Falcon figures; media volumes are working
# volumes for routine culture, not maximum fill.
_VESSELS = {
    '96-well': {'area_cm2': 0.32, 'media_mL': 0.2, 'wells': 96},
    '48-well': {'area_cm2': 0.95, 'media_mL': 0.5, 'wells': 48},
    '24-well': {'area_cm2': 1.9, 'media_mL': 0.5, 'wells': 24},
    '12-well': {'area_cm2': 3.8, 'media_mL': 1.0, 'wells': 12},
    '6-well': {'area_cm2': 9.6, 'media_mL': 2.0, 'wells': 6},
    '10-cm': {'area_cm2': 56.7, 'media_mL': 10.0, 'wells': 1},
    '15-cm': {'area_cm2': 152.0, 'media_mL': 25.0, 'wells': 1},
}

# Spellings people actually type, mapped onto the canonical names above.
_ALIASES = {}
for _name in _VESSELS:
    _stem = _name.split('-')[0]
    _ALIASES[_name] = _name
    _ALIASES[_stem] = _name
    if _name.endswith('-well'):
        _ALIASES[f'{_stem} well'] = _name
        _ALIASES[f'{_stem}well'] = _name
        _ALIASES[f'{_stem}-wells'] = _name
    else:
        _ALIASES[f'{_stem}cm'] = _name
        _ALIASES[f'{_stem} cm'] = _name
        _ALIASES[f'{_stem}cm dish'] = _name
        _ALIASES[f'{_stem} cm dish'] = _name
        _ALIASES[f'{_stem}-cm dish'] = _name
del _name, _stem

# HEK293T seeded one day before transfection, derived from the 96-well number
# this lab already uses (15,000 per well into 0.32 cm^2) so that every other
# format is the same density rather than an independent guess.
SEED_PER_CM2 = 15000 / 0.32

# Format each calculator in edms.bio.transfect was written for
PE3_ANCHOR = '96-well'
# Key used in place of a plate map that would not match the chosen format
LAYOUT_NOTE = 'Plate layout'
VIRUS_ANCHOR = '6-well'

def _round_cells(cells: float) -> int:
    '''
    _round_cells(): round a seeding number to something worth pipetting.

    Parameters:
    cells (float): exact cells per well
    '''
    if cells >= 1e5:  # two significant figures: 2,657,812 -> 2,700,000
        from math import floor, log10
        step = 10 ** (floor(log10(cells)) - 1)
        return int(round(cells / step) * step)
    if cells >= 1e4:
        return int(round(cells / 1000) * 1000)
    return int(round(cells / 100) * 100)

# Formats
def canonical(plate_format: str) -> str:
    '''
    canonical(): resolve a format spelling to its canonical name.

    Parameters:
    plate_format (str): e.g. '96-well', '96', '10cm', '10 cm dish'
    '''
    if not isinstance(plate_format, str):
        raise TypeError(f'plate_format must be a string, got '
                        f'{type(plate_format).__name__}')
    key = ' '.join(plate_format.strip().lower().split())
    if key in _ALIASES:
        return _ALIASES[key]
    raise KeyError(f"Unknown plate format '{plate_format}'. "
                   f"Known formats: {', '.join(names())}")

def names() -> list[str]:
    '''
    names(): canonical format names, smallest growth area first.
    '''
    return sorted(_VESSELS, key=lambda n: _VESSELS[n]['area_cm2'])

def spec(plate_format: str) -> dict:
    '''
    spec(): geometry and seeding recommendation for one format.

    Parameters:
    plate_format (str): plate format name or alias

    Dependencies: canonical()
    '''
    name = canonical(plate_format)
    v = _VESSELS[name]
    return {'plate_format': name,
            'wells': v['wells'],
            'area_cm2': v['area_cm2'],
            'media_mL': v['media_mL'],
            'seed_HEK293T': _round_cells(v['area_cm2'] * SEED_PER_CM2)}

def scale(to: str, frm: str) -> float:
    '''
    scale(): growth-area ratio between two formats.

    Parameters:
    to (str): format being scaled to
    frm (str): format the amount was written for

    Dependencies: spec()
    '''
    return spec(to)['area_cm2'] / spec(frm)['area_cm2']

def table(file: str=None) -> pd.DataFrame:
    '''
    table(): every format as a DataFrame, smallest first.

    Parameters:
    file (str, optional): output file path

    Dependencies: pandas,spec()
    '''
    df = pd.DataFrame([spec(n) for n in names()])
    if file is not None:
        from ..gen import io
        io.save(obj=df, file=file, id=True)
    return df

# Planning
def geometry(plate_format: str, reps: int=1, conditions: int=None,
             file: str=None) -> pd.DataFrame:
    '''
    geometry(): one-row plan of the chosen format, for a protocol step.

    Parameters:
    plate_format (str): plate format name or alias
    reps (int, optional): replicates per condition (Default: 1)
    conditions (int, optional): conditions, to total up wells and media
    file (str, optional): output file path

    Dependencies: pandas,spec()
    '''
    s = spec(plate_format)
    row = {'Format': s['plate_format'],
           'Growth area (cm2)': s['area_cm2'],
           'Media per well (mL)': s['media_mL'],
           'Seed HEK293T per well': s['seed_HEK293T'],
           'Replicates': reps}
    if conditions is not None:
        wells = int(conditions) * int(reps)
        row['Conditions'] = int(conditions)
        row['Wells needed'] = wells
        row['Plates needed'] = round(wells / s['wells'], 2)
        row['Media needed (mL)'] = round(wells * s['media_mL'], 1)
        row['Cells needed'] = wells * s['seed_HEK293T']
    df = pd.DataFrame([row])
    if file is not None:
        from ..gen import io
        io.save(obj=df, file=file, id=True)
    return df

def PE3(plasmids, epegRNAs, ngRNAs, plate_format: str=PE3_ANCHOR,
        epegRNA_ng: float=66, ngRNA_ng: float=22, PE_ng: float=200,
        well_uL: float=10, **kwargs) -> dict:
    '''
    PE3(): edms.bio.transfect:PE3 with DNA and volume rescaled to a format.

    The ng and uL arguments are the per-well amounts for the 96-well plate PE3
    was written for; each is multiplied by the growth-area ratio of the chosen
    format. In a 96-well plate the amounts pass through unchanged.

    Parameters:
    plasmids (DataFrame): nanodrop concentrations for every plasmid
    epegRNAs (DataFrame): epegRNAs from PrimeDesign
    ngRNAs (DataFrame): ngRNAs from PrimeDesign
    plate_format (str, optional): plate format name or alias (Default: 96-well)
    epegRNA_ng (float, optional): epegRNA ng per well at 96-well (Default: 66)
    ngRNA_ng (float, optional): ngRNA ng per well at 96-well (Default: 22)
    PE_ng (float, optional): PE plasmid ng per well at 96-well (Default: 200)
    well_uL (float, optional): transfection mix uL per well at 96-well (Default: 10)
    **kwargs: passed through to edms.bio.transfect:PE3

    The 96-well plate map PE3 returns is dropped for any other format, since it
    describes an 8x12 grid; a note replaces it.

    Dependencies: edms.bio.transfect,scale()
    '''
    name = canonical(plate_format)
    f = scale(name, PE3_ANCHOR)
    tables = transfect.PE3(plasmids=plasmids, epegRNAs=epegRNAs, ngRNAs=ngRNAs,
                           epegRNA_ng=round(epegRNA_ng * f, 2),
                           ngRNA_ng=round(ngRNA_ng * f, 2),
                           PE_ng=round(PE_ng * f, 2),
                           well_uL=round(well_uL * f, 2), **kwargs)

    # transfect.PE3 lays conditions out on a fixed 8x12 grid with the outer
    # wells left empty. That map describes a 96-well plate and nothing else, so
    # for any other format replace it with a note rather than print a layout
    # that does not match the plate on the bench.
    if name != PE3_ANCHOR:
        for key in [k for k in tables if 'plate' in k.lower()]:
            del tables[key]
        tables[f'{LAYOUT_NOTE} ({name})'] = pd.DataFrame([{
            'Note': f'No plate layout: edms.bio.transfect:PE3 lays out an 8x12 '
                    f'grid, which does not describe a {name} plate. The volumes '
                    f'above are scaled correctly; assign wells yourself.'}])
    return tables

def virus(plasmids, plate_format: str=VIRUS_ANCHOR, VSVG_ng: float=750,
          GagPol_ng: float=1500, transfer_ng: float=750, well_uL: float=500,
          P3000_uL_per_well: float=6, L3000_uL_per_well: float=7,
          **kwargs) -> pd.DataFrame:
    '''
    virus(): edms.bio.transfect:virus with DNA and volume rescaled to a format.

    The ng and uL arguments are the per-well amounts for the 6-well plate
    virus() was written for; each is multiplied by the growth-area ratio of the
    chosen format. In a 6-well plate the amounts pass through unchanged.

    Parameters:
    plasmids (DataFrame): nanodrop concentrations for every plasmid
    plate_format (str, optional): plate format name or alias (Default: 6-well)
    VSVG_ng (float, optional): VSVG ng per well at 6-well (Default: 750)
    GagPol_ng (float, optional): GagPol ng per well at 6-well (Default: 1500)
    transfer_ng (float, optional): transfer plasmid ng per well at 6-well (Default: 750)
    well_uL (float, optional): transfection mix uL per well at 6-well (Default: 500)
    P3000_uL_per_well (float, optional): P3000 uL per well at 6-well (Default: 6)
    L3000_uL_per_well (float, optional): L3000 uL per well at 6-well (Default: 7)
    **kwargs: passed through to edms.bio.transfect:virus

    Dependencies: edms.bio.transfect,scale()
    '''
    f = scale(plate_format, VIRUS_ANCHOR)
    if hasattr(plasmids, 'copy'):  # virus() writes columns onto the frame
        plasmids = plasmids.copy()
    return transfect.virus(plasmids=plasmids,
                           VSVG_ng=round(VSVG_ng * f, 2),
                           GagPol_ng=round(GagPol_ng * f, 2),
                           transfer_ng=round(transfer_ng * f, 2),
                           well_uL=round(well_uL * f, 2),
                           P3000_uL_per_well=round(P3000_uL_per_well * f, 2),
                           L3000_uL_per_well=round(L3000_uL_per_well * f, 2),
                           **kwargs)
