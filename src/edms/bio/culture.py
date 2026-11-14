'''
Module: culture.py
Author: Marc Zepeda
Created: 2026-11-13
Description: Cell culture master mixes

The same arithmetic turns up at every step between transduction and harvest:
a compound at a stock concentration has to end up at a desired concentration
in a well, after cells in their own volume of media are added on top. Polybrene
at transduction, puromycin at selection, DHT or ARV-110 at treatment are all
the same calculation with different units.

Because the drug goes into only part of the final volume, it has to be made up
over-concentrated by the ratio of the final volume to the drug volume. Getting
that ratio wrong is the quiet failure mode: the cells see less drug than the
protocol claims, and an incomplete selection looks like a biological result.

Usage:
[Master mixes]
- drug_mm(): master mix for a compound added to wells ahead of cells
- cells_mm(): volume to spin down and resuspend to seed wells at a density
- transduction(): the pair of mixes a transduction needs

[Seeding from confluency]
- confluent_cells(): cells a vessel holds at a given confluency
- per_mL(): estimate a harvest's concentration from how confluent it looked
- seed_plan(): spin volume per sample for one arm of an experiment
- competitive_growth(): both arms of a transduction + competitive growth setup

[Clone isolation]
- sort_plan(): plates to sort and conditioned media to prepare for a clone hunt
- freeze_down(): vials and freezing medium from a harvest
'''
import math
import pandas as pd

from ..gen import io

# Cells in a 100% confluent vessel, from Thermo Fisher's 'Useful Numbers for
# Cell Culture' table. Measured with HeLa, so a bigger or smaller line will
# sit above or below these; they are a starting point, not a property of the
# plastic.
CONFLUENT_CELLS = {
    '96-well': 0.04e6,
    '48-well': 0.12e6,
    '24-well': 0.24e6,
    '12-well': 0.5e6,
    '6-well': 1.2e6,
    '6-cm': 3.2e6,
    '10-cm': 8.8e6,
    '15-cm': 20.0e6,
}

# Assumed concentration (cells/mL) of a harvest from a vessel that was 100%
# confluent, trypsinised and resuspended the usual way. This conflates two
# things -- how many cells the vessel held, and how much media they went into
# -- which is why the same entry can give two different values for one vessel:
# MUZ361 used 1.5e6 for a 15 cm dish in one week and 2.0e6 the next.
#
# Both of those are the same cell count at a different resuspension volume.
# The values below are 10/30/26's, and they are exactly the Thermo confluent
# counts above at a 1 mL resuspension for a 6-well and a 10 mL resuspension
# for a 15 cm dish. Pass harvest_mL to per_mL() and it works the concentration
# out that way instead, for any vessel.
CONFLUENT_PER_ML = {
    '6-well': CONFLUENT_CELLS['6-well'] / 1.0,
    '15-cm': CONFLUENT_CELLS['15-cm'] / 10.0,
}

def drug_mm(stock: float, desired: float, wells: int, drug_media_mL: float=1.0,
            cell_media_mL: float=1.0, mm_x: float=1.1, drug: str='Drug',
            unit: str='ug/mL', file: str=None) -> pd.DataFrame:
    '''
    drug_mm(): master mix for a compound added to wells ahead of cells.

    The drug is made up in drug_media_mL per well and the cells arrive in
    cell_media_mL on top, so the mix is concentrated by
    (drug_media_mL + cell_media_mL) / drug_media_mL to land at `desired` in the
    final volume.

    Parameters:
    stock (float): stock concentration, in `unit`
    desired (float): concentration wanted in the well after cells are added
    wells (int): wells to make up
    drug_media_mL (float, optional): drug media per well in mL (Default: 1.0)
    cell_media_mL (float, optional): media the cells arrive in, mL (Default: 1.0)
    mm_x (float, optional): master mix multiplier (Default: 1.1)
    drug (str, optional): name used in the table (Default: 'Drug')
    unit (str, optional): concentration unit for stock and desired (Default: 'ug/mL')
    file (str, optional): output file path

    Dependencies: pandas,io
    '''
    if drug_media_mL <= 0:
        raise ValueError(f'drug_media_mL must be greater than 0, got {drug_media_mL}')
    if stock <= 0:
        raise ValueError(f'stock must be greater than 0, got {stock}')
    if desired > stock:
        raise ValueError(f'desired ({desired} {unit}) is above the stock '
                         f'concentration ({stock} {unit}); it cannot be reached '
                         f'by dilution')

    final_mL = drug_media_mL + cell_media_mL
    conc_mm = desired * final_mL / drug_media_mL
    total_mL = drug_media_mL * wells * mm_x
    stock_mL = total_mL * conc_mm / stock
    media_mL = total_mL - stock_mL

    rows = [(f'{drug}, stock', stock, unit),
            (f'{drug}, desired', desired, unit),
            ('Drug media / well', drug_media_mL, 'mL'),
            ('Cell media / well', cell_media_mL, 'mL'),
            ('Wells', wells, 'wells'),
            ('MM', mm_x, 'x'),
            (f'{drug}, MM', round(conc_mm, 4), unit),
            (f'{drug} stock, MM', round(stock_mL, 4), 'mL'),
            ('Media, MM', round(media_mL, 4), 'mL')]
    df_drug = pd.DataFrame({'Component': [r[0] for r in rows],
                            'Value': pd.Series([r[1] for r in rows], dtype=object),
                            'Unit': [r[2] for r in rows]})
    if file is not None:
        io.save(obj=df_drug, file=file, id=True)
    return df_drug

def cells_mm(flow_per_mL: float, cells_per_well: float, wells: int,
             cell_media_mL: float=1.0, mm_x: float=1.1,
             file: str=None) -> pd.DataFrame:
    '''
    cells_mm(): volume to spin down and resuspend to seed wells at a density.

    Parameters:
    flow_per_mL (float): measured concentration of the culture, cells/mL
    cells_per_well (float): cells wanted per well
    wells (int): wells to seed
    cell_media_mL (float, optional): volume the cells go into per well, mL (Default: 1.0)
    mm_x (float, optional): master mix multiplier (Default: 1.1)
    file (str, optional): output file path

    Dependencies: pandas,io
    '''
    if flow_per_mL <= 0:
        raise ValueError(f'flow_per_mL must be greater than 0, got {flow_per_mL}')

    cells_needed = cells_per_well * wells * mm_x
    spin_mL = cells_needed / flow_per_mL
    resuspend_mL = cell_media_mL * wells * mm_x

    rows = [('Cells, flow', flow_per_mL, 'cells/mL'),
            ('Cells, desired', cells_per_well, 'cells/well'),
            ('Cell media / well', cell_media_mL, 'mL'),
            ('Wells', wells, 'wells'),
            ('MM', mm_x, 'x'),
            ('Cells needed, MM', round(cells_needed, 1), 'cells'),
            ('Resuspension, MM', round(cells_needed / resuspend_mL, 1), 'cells/mL'),
            ('Culture to spin, MM', round(spin_mL, 4), 'mL'),
            ('Resuspend in, MM', round(resuspend_mL, 4), 'mL')]
    df_cells = pd.DataFrame({'Component': [r[0] for r in rows],
                             'Value': pd.Series([r[1] for r in rows], dtype=object),
                             'Unit': [r[2] for r in rows]})
    if file is not None:
        io.save(obj=df_cells, file=file, id=True)
    return df_cells

def transduction(wells: int, flow_per_mL: float, cells_per_well: float,
                 virus_mL_per_well: float=1.0, cell_media_mL: float=1.0,
                 polybrene_stock: float=10000, polybrene_desired: float=8,
                 polybrene_unit: str='ug/mL', mm_x_polybrene: float=1.2,
                 mm_x_cells: float=1.1, file: str=None) -> dict[str, pd.DataFrame]:
    '''
    transduction(): the pair of mixes a transduction needs.

    Virus media goes onto the plate first, then cells resuspended in polybrene
    media on top, so polybrene is made up over the virus volume as well as its
    own. The two mixes carry different multipliers because the polybrene media
    is cheap to overmake and the cells are not.

    Parameters:
    wells (int): wells to transduce
    flow_per_mL (float): measured concentration of the culture, cells/mL
    cells_per_well (float): cells wanted per well
    virus_mL_per_well (float, optional): virus media per well, mL (Default: 1.0)
    cell_media_mL (float, optional): polybrene media the cells arrive in, mL (Default: 1.0)
    polybrene_stock (float, optional): polybrene stock (Default: 10000)
    polybrene_desired (float, optional): polybrene wanted in the well (Default: 8)
    polybrene_unit (str, optional): polybrene concentration unit (Default: 'ug/mL')
    mm_x_polybrene (float, optional): polybrene master mix multiplier (Default: 1.2)
    mm_x_cells (float, optional): cell master mix multiplier (Default: 1.1)
    file (str, optional): output file path

    Dependencies: pandas,io,drug_mm(),cells_mm()
    '''
    # Polybrene is made up in the cell media but has to cover the virus volume
    # too, so the virus volume is what the cells are "added on top of" here.
    polybrene = drug_mm(stock=polybrene_stock, desired=polybrene_desired,
                        wells=wells, drug_media_mL=cell_media_mL,
                        cell_media_mL=virus_mL_per_well, mm_x=mm_x_polybrene,
                        drug='Polybrene', unit=polybrene_unit)
    cells = cells_mm(flow_per_mL=flow_per_mL, cells_per_well=cells_per_well,
                     wells=wells, cell_media_mL=cell_media_mL, mm_x=mm_x_cells)
    tables = {'Polybrene': polybrene, 'Cells': cells}
    if file is not None:
        io.save(obj=tables, file=file, id=True)
    return tables

# Seeding from confluency
def confluent_cells(vessel: str, confluency: float=100, recovery: float=1.0,
                    confluent_cells_: dict=None) -> float:
    '''
    confluent_cells(): cells a vessel holds at a given confluency.

    Unlike per_mL() this does not depend on what the harvest is resuspended
    in, so it is the right thing to ask when the question is how many cells
    you have rather than how concentrated they are.

    The counts are Thermo's HeLa figures. A line that is not HeLa will not
    match: MUZ361 froze 12 vials of a million cells from a 15 cm dish at 80%,
    where these numbers predict 16, so LNCaP comes in around three quarters of
    a HeLa dish once trypsinisation and spin losses are taken off. That is
    what `recovery` is for.

    Parameters:
    vessel (str): culture vessel, e.g. '6-well', '15cm', '150mm'
    confluency (float, optional): percent confluent, 0-100 (Default: 100)
    recovery (float, optional): fraction of the predicted cells actually
        recovered, above 0 and at most 1 (Default: 1.0)
    confluent_cells_ (dict, optional): confluent counts per vessel, overriding
        CONFLUENT_CELLS for the vessels it names

    Dependencies: edms.work.plates,CONFLUENT_CELLS
    '''
    from ..work import plates  # local: edms.work.plates imports edms.bio

    if not 0 < confluency <= 100:
        raise ValueError(f'confluency must be a percent above 0 and at most '
                         f'100, got {confluency}')
    if not 0 < recovery <= 1:
        raise ValueError(f'recovery is a fraction above 0 and at most 1, got '
                         f'{recovery}')
    name = plates.canonical(vessel)
    table = dict(CONFLUENT_CELLS)
    for key, value in (confluent_cells_ or {}).items():
        table[plates.canonical(key)] = value
    if name not in table:
        known = ', '.join(sorted(table)) or '(none)'
        raise KeyError(f"No confluent cell count for '{name}'. Known vessels: "
                       f'{known}. Pass confluent_cells_ to add it.')
    return table[name] * confluency / 100 * recovery

def per_mL(vessel: str, confluency: float, harvest_mL: float=None,
           confluent_per_mL: dict=None) -> float:
    '''
    per_mL(): estimate a harvest's concentration from how confluent it looked.

    Counting every flask is the accurate way to do this and the slow way. When
    fifteen clones come off the incubator at once, the usual shortcut is to
    read confluency off the microscope and scale a known full-plate yield. The
    estimate is linear in confluency, which is close enough over the 40-100%
    range people actually seed from.

    Parameters:
    vessel (str): culture vessel the harvest came from, e.g. '6-well', '15cm'
    confluency (float): percent confluent, 0-100
    harvest_mL (float, optional): volume the harvest is resuspended in. Given
        this, the concentration is worked out from the vessel's confluent cell
        count, which covers every vessel rather than only the calibrated ones
        (Default: None)
    confluent_per_mL (dict, optional): cells/mL at 100% confluency per vessel,
        overriding CONFLUENT_PER_ML for the vessels it names

    Dependencies: edms.work.plates,confluent_cells(),CONFLUENT_PER_ML
    '''
    from ..work import plates  # local: edms.work.plates imports edms.bio

    if not 0 < confluency <= 100:
        raise ValueError(f'confluency must be a percent above 0 and at most '
                         f'100, got {confluency}')
    if harvest_mL is not None:
        if harvest_mL <= 0:
            raise ValueError(f'harvest_mL must be greater than 0, got {harvest_mL}')
        return confluent_cells(vessel=vessel, confluency=confluency) / harvest_mL
    name = plates.canonical(vessel)
    table = dict(CONFLUENT_PER_ML)
    for key, value in (confluent_per_mL or {}).items():
        table[plates.canonical(key)] = value
    if name not in table:
        known = ', '.join(sorted(table)) or '(none)'
        raise KeyError(f"No confluent density for '{name}'. Known vessels: "
                       f'{known}. Pass harvest_mL and it is worked out from '
                       f'the confluent cell count instead, or pass '
                       f'confluent_per_mL, or give the measured concentration '
                       f'per sample.')
    return table[name] * confluency / 100

def seed_plan(df: pd.DataFrame, cells_per_well: float, wells: int,
              cell_media_mL: float=0.25, mm_x: float=1.2, arm: str='Seed',
              confluency_col: str='Confluency %', vessel_col: str='Vessel',
              per_mL_col: str='Cells, flow', scale_col: str=None,
              confluent_per_mL: dict=None, file: str=None) -> pd.DataFrame:
    '''
    seed_plan(): spin volume per sample for one arm of an experiment.

    One row per sample, each at its own confluency in its own vessel, each
    going into the same number of wells. Every sample needs the same number of
    cells, so the only thing that varies is how much of its harvest to spin
    down. A sample whose concentration was actually measured can carry it in
    per_mL_col and that number is used instead of the confluency estimate.

    Called more than once on the same table with different `arm` names, the
    columns accumulate and a 'Sum spin, mL' column appears: the total to take
    from each sample when both arms are set up off one harvest.

    Parameters:
    df (dataframe): one row per sample, with confluency and vessel columns
    cells_per_well (float): cells wanted per well
    wells (int): wells this arm seeds per sample
    cell_media_mL (float, optional): volume the cells go into per well, mL (Default: 0.25)
    mm_x (float, optional): master mix multiplier (Default: 1.2)
    arm (str, optional): name for this arm's columns (Default: 'Seed')
    confluency_col (str, optional): percent confluent column (Default: 'Confluency %')
    vessel_col (str, optional): vessel column (Default: 'Vessel')
    per_mL_col (str, optional): optional measured cells/mL column, used in
        place of the estimate where it holds a number (Default: 'Cells, flow')
    scale_col (str, optional): optional column of per-sample multipliers on the
        cells needed, for a sample deliberately seeded heavier. It applies to
        whichever arm this call is planning, so a sample seeded heavier in one
        arm only means calling seed_plan() per arm rather than through
        competitive_growth() (Default: None)
    confluent_per_mL (dict, optional): cells/mL at 100% confluency per vessel
    file (str, optional): output file path

    Dependencies: pandas,io,per_mL()
    '''
    for col in (confluency_col, vessel_col):
        if col not in df.columns:
            raise KeyError(f"seed_plan() needs a '{col}' column; df has "
                           f"{', '.join(map(str, df.columns))}")
    if wells <= 0:
        raise ValueError(f'wells must be greater than 0, got {wells}')

    df_plan = df.copy()
    cells_needed = cells_per_well * wells * mm_x
    resuspend_mL = cell_media_mL * wells * mm_x

    spins, resuspends = [], []
    for _, row in df_plan.iterrows():
        measured = row.get(per_mL_col) if per_mL_col in df_plan.columns else None
        if pd.notna(measured) and measured not in ('', None):
            conc = float(measured)
            if conc <= 0:
                raise ValueError(f"'{per_mL_col}' must be greater than 0, got {conc}")
        else:
            conc = per_mL(vessel=row[vessel_col], confluency=float(row[confluency_col]),
                          confluent_per_mL=confluent_per_mL)
        scale = float(row[scale_col]) if scale_col and pd.notna(row.get(scale_col)) else 1.0
        spins.append(round(cells_needed * scale / conc, 4))
        resuspends.append(round(resuspend_mL * scale, 4))

    df_plan[f'{arm} spin, mL'] = spins
    df_plan[f'{arm} resuspend, mL'] = resuspends

    spin_cols = [c for c in df_plan.columns if c.endswith(' spin, mL')]
    if len(spin_cols) > 1:
        df_plan['Sum spin, mL'] = df_plan[spin_cols].sum(axis=1).round(4)
        # Keep the total last however many arms have been added
        df_plan = df_plan[[c for c in df_plan.columns if c != 'Sum spin, mL']
                          + ['Sum spin, mL']]
    if file is not None:
        io.save(obj=df_plan, file=file, id=True)
    return df_plan

def competitive_growth(df: pd.DataFrame, cg_cells_per_well: float=1.2e5,
                       cg_wells: int=3, transduce_cells_per_well: float=1.8e5,
                       transduce_wells: int=3, cell_media_mL: float=0.25,
                       mm_x: float=1.2, competitor: str='LCG',
                       competitor_vessel: str='15-cm',
                       competitor_confluency: float=100,
                       competitor_per_mL: float=None, competitor_mm_x: float=1.2,
                       extra_samples: int=1, confluency_col: str='Confluency %',
                       vessel_col: str='Vessel', per_mL_col: str='Cells, flow',
                       scale_col: str=None, confluent_per_mL: dict=None,
                       file: str=None) -> dict[str, pd.DataFrame]:
    '''
    competitive_growth(): both arms of a transduction + competitive growth setup.

    One harvest per sample feeds two plates: a transduction arm, and a
    competitive growth arm seeded 50:50 against a common competitor line. The
    arms want different numbers of cells, so the table gives each arm's spin
    volume and the total to take from each sample.

    The competitor is prepared once for every sample at a multiplier of its
    own, plus extra_samples worth for the day 0 mixing readout, since getting
    it wrong wastes the whole comparison rather than one well.

    Set transduce_wells to 0 for a competitive growth arm on its own.

    Parameters:
    df (dataframe): one row per sample, with confluency and vessel columns
    cg_cells_per_well (float, optional): competitive growth cells per well from
        the sample; the competitor matches it (Default: 1.2e5)
    cg_wells (int, optional): competitive growth wells per sample (Default: 3)
    transduce_cells_per_well (float, optional): transduction cells per well (Default: 1.8e5)
    transduce_wells (int, optional): transduction wells per sample, 0 to skip (Default: 3)
    cell_media_mL (float, optional): volume the cells go into per well, mL (Default: 0.25)
    mm_x (float, optional): master mix multiplier (Default: 1.2)
    competitor (str, optional): competitor line name (Default: 'LCG')
    competitor_vessel (str, optional): vessel the competitor came from (Default: '15-cm')
    competitor_confluency (float, optional): percent confluent (Default: 100)
    competitor_per_mL (float, optional): measured cells/mL, used in place of
        the confluency estimate (Default: None)
    competitor_mm_x (float, optional): extra multiplier on the competitor (Default: 1.2)
    extra_samples (int, optional): samples' worth of competitor beyond the
        table, for the day 0 mix (Default: 1)
    confluency_col (str, optional): percent confluent column (Default: 'Confluency %')
    vessel_col (str, optional): vessel column (Default: 'Vessel')
    per_mL_col (str, optional): optional measured cells/mL column (Default: 'Cells, flow')
    scale_col (str, optional): optional per-sample multiplier column (Default: None)
    confluent_per_mL (dict, optional): cells/mL at 100% confluency per vessel
    file (str, optional): output file path

    Dependencies: pandas,io,per_mL(),seed_plan()
    '''
    if len(df) == 0:
        raise ValueError('competitive_growth() needs at least one sample')

    shared = dict(cell_media_mL=cell_media_mL, mm_x=mm_x,
                  confluency_col=confluency_col, vessel_col=vessel_col,
                  per_mL_col=per_mL_col, scale_col=scale_col,
                  confluent_per_mL=confluent_per_mL)
    df_plan = df
    if transduce_wells:
        df_plan = seed_plan(df=df_plan, cells_per_well=transduce_cells_per_well,
                            wells=transduce_wells, arm='Transduce', **shared)
    df_plan = seed_plan(df=df_plan, cells_per_well=cg_cells_per_well,
                        wells=cg_wells, arm='CG', **shared)

    # The competitor is made up once for every sample, plus the day 0 mix
    samples = len(df) + extra_samples
    cells_needed = cg_cells_per_well * cg_wells * mm_x * samples * competitor_mm_x
    resuspend_mL = cell_media_mL * cg_wells * mm_x * samples * competitor_mm_x
    if competitor_per_mL is not None:
        if competitor_per_mL <= 0:
            raise ValueError(f'competitor_per_mL must be greater than 0, got '
                             f'{competitor_per_mL}')
        conc = competitor_per_mL
        source = 'measured'
    else:
        conc = per_mL(vessel=competitor_vessel, confluency=competitor_confluency,
                      confluent_per_mL=confluent_per_mL)
        source = f'{competitor_confluency:g}% of a {competitor_vessel}'
    rows = [(f'{competitor}, source', source, '-'),
            (f'{competitor}, cells', round(conc, 1), 'cells/mL'),
            ('Cells, desired', cg_cells_per_well, 'cells/well'),
            ('Wells / sample', cg_wells, 'wells'),
            ('Samples', samples, 'samples'),
            ('MM', round(mm_x * competitor_mm_x, 4), 'x'),
            (f'{competitor} needed, MM', round(cells_needed, 1), 'cells'),
            (f'{competitor} to spin, MM', round(cells_needed / conc, 4), 'mL'),
            ('Resuspend in, MM', round(resuspend_mL, 4), 'mL')]
    df_competitor = pd.DataFrame({'Component': [r[0] for r in rows],
                                  'Value': pd.Series([r[1] for r in rows], dtype=object),
                                  'Unit': [r[2] for r in rows]})

    tables = {'Seed plan': df_plan, competitor: df_competitor}
    if file is not None:
        io.save(obj=tables, file=file, id=True)
    return tables

# Clone isolation
def sort_plan(target_clones: int, plating_efficiency: float=0.1,
              plate_format: str='96-well', conditioned_fraction: float=0.5,
              media_mL: float=None, file: str=None) -> pd.DataFrame:
    '''
    sort_plan(): plates to sort and conditioned media to prepare for a clone hunt.

    Single cells sorted one per well do not all become colonies, so the number
    of plates to sort is the number of clones wanted divided by the fraction
    that survive. That fraction is the whole calculation and it is
    line-specific: LNCaP is at the unforgiving end of it.

    No entry here records a measured plating efficiency, so the default is a
    placeholder rather than a figure to trust. Sort one plate, count the wells
    that grew, and put your own number in before planning the real sort.

    The conditioned media is the part that has to be ready in advance. A sort
    that happens before it is made either waits or goes out into plain media,
    and clones that will not start in plain media are lost at that point
    rather than later.

    Parameters:
    target_clones (int): clones wanted at the end
    plating_efficiency (float, optional): fraction of sorted wells that become
        a colony, between 0 and 1 (Default: 0.1)
    plate_format (str, optional): format sorted into (Default: '96-well')
    conditioned_fraction (float, optional): fraction of the well volume that is
        conditioned media (Default: 0.5)
    media_mL (float, optional): media per well, mL; taken from the format when
        left unset (Default: None)
    file (str, optional): output file path

    Dependencies: math,pandas,io,edms.work.plates
    '''
    from ..work import plates  # local: edms.work.plates imports edms.bio

    if target_clones < 1:
        raise ValueError(f'target_clones must be at least 1, got {target_clones}')
    if not 0 < plating_efficiency <= 1:
        raise ValueError(f'plating_efficiency is a fraction above 0 and at '
                         f'most 1, got {plating_efficiency}')
    if not 0 <= conditioned_fraction <= 1:
        raise ValueError(f'conditioned_fraction is a fraction between 0 and 1, '
                         f'got {conditioned_fraction}')

    s = plates.spec(plate_format)
    per_well_mL = s['media_mL'] if media_mL is None else media_mL
    wells_needed = math.ceil(target_clones / plating_efficiency)
    plate_count = math.ceil(wells_needed / s['wells'])
    wells_sorted = plate_count * s['wells']
    total_mL = wells_sorted * per_well_mL

    rows = [('Clones wanted', target_clones, 'clones'),
            ('Plating efficiency', plating_efficiency, 'fraction'),
            ('Format', s['plate_format'], '-'),
            ('Wells needed', wells_needed, 'wells'),
            ('Plates to sort', plate_count, 'plates'),
            ('Wells sorted', wells_sorted, 'wells'),
            ('Clones expected', round(wells_sorted * plating_efficiency, 1), 'clones'),
            ('Media per well', per_well_mL, 'mL'),
            ('Media total', round(total_mL, 1), 'mL'),
            ('Conditioned media', round(total_mL * conditioned_fraction, 1), 'mL'),
            ('Fresh media', round(total_mL * (1 - conditioned_fraction), 1), 'mL')]
    df_sort = pd.DataFrame({'Component': [r[0] for r in rows],
                            'Value': pd.Series([r[1] for r in rows], dtype=object),
                            'Unit': [r[2] for r in rows]})
    if file is not None:
        io.save(obj=df_sort, file=file, id=True)
    return df_sort

def freeze_down(harvest_mL: float, cells_per_mL: float=None, vessel: str=None,
                confluency: float=None, recovery: float=1.0,
                cells_per_vial: float=1e6, vial_mL: float=1.0,
                dmso_pct: float=10, confluent_cells_: dict=None,
                file: str=None) -> pd.DataFrame:
    '''
    freeze_down(): vials and freezing medium from a harvest.

    Give it a measured concentration, or a vessel and a confluency to estimate
    one. Vials round down, because a part-filled vial is not a vial you would
    thaw.

    How many cells a dish holds does not depend on what you resuspend them in,
    so the estimate goes through confluent_cells() and harvest_mL only sets the
    concentration you end up working at.

    MUZ361 froze 12 vials from a 15 cm dish at 75-80% confluent, twice, where
    Thermo's HeLa counts predict 16. That is what recovery is for: 0.75 gets
    those entries back, and is a reasonable starting point for LNCaP once
    trypsinisation and spin losses are taken off.

    Parameters:
    harvest_mL (float): volume the pellet is resuspended in before counting, mL
    cells_per_mL (float, optional): measured concentration (Default: None)
    vessel (str, optional): vessel the harvest came from, for the estimate
        (Default: None)
    confluency (float, optional): percent confluent, for the estimate (Default: None)
    recovery (float, optional): fraction of the predicted cells actually
        recovered (Default: 1.0)
    cells_per_vial (float, optional): cells to put in each vial (Default: 1e6)
    vial_mL (float, optional): freezing medium per vial, mL (Default: 1.0)
    dmso_pct (float, optional): DMSO in the freezing medium, percent (Default: 10)
    confluent_cells_ (dict, optional): confluent counts per vessel
    file (str, optional): output file path

    Dependencies: pandas,io,confluent_cells()
    '''
    if harvest_mL <= 0:
        raise ValueError(f'harvest_mL must be greater than 0, got {harvest_mL}')
    if cells_per_vial <= 0:
        raise ValueError(f'cells_per_vial must be greater than 0, got {cells_per_vial}')
    if not 0 <= dmso_pct < 100:
        raise ValueError(f'dmso_pct must be at least 0 and below 100, got {dmso_pct}')

    if cells_per_mL is not None:
        if cells_per_mL <= 0:
            raise ValueError(f'cells_per_mL must be greater than 0, got {cells_per_mL}')
        conc, source = float(cells_per_mL), 'measured'
        total = conc * harvest_mL
    elif vessel is not None and confluency is not None:
        from ..work import plates  # local: edms.work.plates imports edms.bio
        total = confluent_cells(vessel=vessel, confluency=confluency,
                                recovery=recovery,
                                confluent_cells_=confluent_cells_)
        conc = total / harvest_mL
        source = f'{confluency:g}% of a {plates.canonical(vessel)}'
        if recovery != 1.0:
            source += f' at {recovery:g} recovery'
    else:
        raise ValueError('freeze_down() needs either cells_per_mL, or both '
                         'vessel and confluency to estimate one')

    vials = int(total // cells_per_vial)
    if vials < 1:
        raise ValueError(
            f'{total:,.0f} cells is less than the {cells_per_vial:,.0f} wanted '
            f'in one vial. Grow the clone further, or lower cells_per_vial.')
    medium_mL = vials * vial_mL
    leftover = total - vials * cells_per_vial

    rows = [('Cells, source', source, '-'),
            ('Cells', round(conc, 1), 'cells/mL'),
            ('Harvest', harvest_mL, 'mL'),
            ('Cells, total', round(total, 1), 'cells'),
            ('Cells / vial', cells_per_vial, 'cells'),
            ('Vials', vials, 'vials'),
            ('Cells left over', round(leftover, 1), 'cells'),
            ('Freezing medium', round(medium_mL, 2), 'mL'),
            ('DMSO', round(medium_mL * dmso_pct / 100, 2), 'mL'),
            ('Media / FBS', round(medium_mL * (1 - dmso_pct / 100), 2), 'mL')]
    df_freeze = pd.DataFrame({'Component': [r[0] for r in rows],
                              'Value': pd.Series([r[1] for r in rows], dtype=object),
                              'Unit': [r[2] for r in rows]})
    if file is not None:
        io.save(obj=df_freeze, file=file, id=True)
    return df_freeze

