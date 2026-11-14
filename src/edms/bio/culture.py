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
'''
import pandas as pd

from ..gen import io

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
