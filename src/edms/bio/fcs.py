''' 
Module: fcs.py
Author: Marc Zepeda
Created: 2026-11-09
Description: Flow Cytometry Standard (FCS) files

Usage:
[NovoCyte .ncf file retrieval]
- ncf_samples(): retrieve sample metadata from a NovoCyte .ncf file as a dataframe
- ncf_events(): retrieve one sample's event matrix from a NovoCyte .ncf file as a dataframe

[FCS file generation]
- fcs_write(): write an event dataframe to an FCS 3.1 file
- ncf_to_fcs(): convert every sample holding data in a NovoCyte .ncf file to FCS 3.1 files

[Flow cytometry plots]
- offset_summary(): draw a half-offset distribution of a flow channel beside its summary statistic
                    and, optionally, that statistic as a % difference from reference sample(s)

Notes:
A NovoCyte .ncf file is a SQLite database. Event data lives in SampleDataData as
float32 little-endian blobs stored column-major within each chunk, with the chunks
concatenated in [Order]. FCS requires the opposite layout (row-major, event-major),
so the matrix is transposed on write.

The NovoCyte database names detectors by laser and position (BL1, YL2) whereas the
NovoExpress FCS export names them by laser letter and emission wavelength (B530,
Y615). NCF_RENAME holds the two mappings confirmed against the NovoCyte 3000RYB
Operator's Guide (Table 20); the original name is kept in $PnS.
'''
# Import packages
import array
import math
import os
import re
import sqlite3
import sys
import numpy as np
import pandas as pd

from ..gen import io
from ..gen import plot as p

# Detector renaming: BL1 -> B530 (488 nm, 530/30; FITC/GFP),
#                    YL2 -> Y615 (561 nm, 615/20; PE-Texas Red/mCherry)
NCF_RENAME = {'BL1': 'B530', 'YL2': 'Y615'}

CHANNEL_RE = re.compile(r'^([A-Za-z]+\d+)(-[AHW])?$')
FCS_DELIMS = '/|\\!#%&*+:;?@^~'
MONTHS = 'JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC'.split()

# NovoCyte .ncf file retrieval
def ncf_samples(pt: str) -> pd.DataFrame:
    ''' 
    ncf_samples(): retrieve sample metadata from a NovoCyte .ncf file as a dataframe

    Parameters:
    pt (str): path to NovoCyte .ncf file

    Dependencies: sqlite3 & pandas
    '''
    if not os.path.isfile(pt):
        raise FileNotFoundError(f'No such .ncf file: {pt}')
    con = sqlite3.connect(f'file:{pt}?mode=ro', uri=True)
    try:
        rows = con.execute('''
            SELECT d.SD_ID, c.SC_ID, c.WellID, c.SampleName, d.Events,
                   d.AcquisitionTime, d.Duration, d.Volume, d.Operator,
                   c.ParameterNames
            FROM SampleData d
            JOIN SampleConfig c ON c.SC_ID = d.SC_ID
            ORDER BY d.SD_ID
        ''').fetchall()
    finally:
        con.close()
    cols = ['SD_ID', 'SC_ID', 'WellID', 'SampleName', 'Events', 'AcquisitionTime',
            'Duration', 'Volume', 'Operator', 'ParameterNames']
    df = pd.DataFrame(rows, columns=cols)
    # Trailing tab in ParameterNames yields a phantom channel if not stripped.
    df['Parameters'] = [(p or '').strip().split('\t') for p in df['ParameterNames']]
    df['PAR'] = [len(p) for p in df['Parameters']]
    return df.drop(columns=['ParameterNames'])

def ncf_events(pt: str, sd_id: int, parameters: list=None, rename: dict | bool = True) -> pd.DataFrame:
    ''' 
    ncf_events(): retrieve one sample's event matrix from a NovoCyte .ncf file as a dataframe

    Parameters:
    pt (str): path to NovoCyte .ncf file
    sd_id (int): SD_ID of the sample to retrieve
    parameters (list, optional): channel names in file order (Default: None; read from SampleConfig)
    rename (dict | bool, optional): detector renaming map, True for NCF_RENAME, or False for none (Default: True)

    Dependencies: array, sqlite3, sys, & pandas
    '''
    mapping = NCF_RENAME if rename is True else ({} if not rename else dict(rename))
    con = sqlite3.connect(f'file:{pt}?mode=ro', uri=True)
    try:
        if parameters is None:
            row = con.execute('''SELECT c.ParameterNames FROM SampleConfig c
                                 JOIN SampleData d ON d.SC_ID = c.SC_ID
                                 WHERE d.SD_ID = ?''', (sd_id,)).fetchone()
            if row is None:
                raise ValueError(f'SD_ID {sd_id} not found in {os.path.basename(pt)}')
            parameters = (row[0] or '').strip().split('\t')
        par = len(parameters)
        if par == 0 or parameters == ['']:
            raise ValueError(f'SD_ID {sd_id}: no parameter names')

        cols = [array.array('f') for _ in range(par)]
        for (blob,) in con.execute('''SELECT Data FROM SampleDataData
                                      WHERE SD_ID = ? ORDER BY [Order]''', (sd_id,)):
            chunk = array.array('f')
            chunk.frombytes(blob)
            if sys.byteorder != 'little':
                chunk.byteswap()
            if len(chunk) % par:
                raise ValueError(f'SD_ID {sd_id}: chunk of {len(chunk)} floats is '
                                 f'not a multiple of {par} parameters')
            n = len(chunk) // par
            for i in range(par): # Column-major within each chunk
                cols[i].extend(chunk[i * n:(i + 1) * n])
    finally:
        con.close()

    names = [_rename_channel(p, mapping) for p in parameters]
    df = pd.DataFrame({name: pd.Series(col, dtype='float32') for name, col in zip(names, cols)})
    df.attrs['originals'] = {new: old for new, old in zip(names, parameters) if new != old}
    return df

def _rename_channel(name: str, mapping: dict) -> str:
    ''' 
    _rename_channel(): apply a detector renaming map to one channel name

    Parameters:
    name (str): channel name (e.g., BL1-A)
    mapping (dict): detector renaming map (e.g., NCF_RENAME)

    Dependencies: re
    '''
    match = CHANNEL_RE.match(name)
    if match and match.group(1) in mapping:
        return mapping[match.group(1)] + (match.group(2) or '')
    return name

# FCS file generation
def fcs_write(df: pd.DataFrame, file: str, keywords: dict=None, originals: dict=None):
    ''' 
    fcs_write(): write an event dataframe to an FCS 3.1 file

    Parameters:
    df (dataframe): events (rows) by channels (columns); every column must be numeric
    file (str): output file path (.fcs)
    keywords (dict, optional): extra FCS keywords, $-prefixed for standard ones (Default: None)
    originals (dict, optional): channel name -> original detector name, stored in $PnS (Default: None)

    Dependencies: array, math, os, sys, & pandas
    '''
    names = list(df.columns)
    par, tot = len(names), len(df)
    if par == 0 or tot == 0:
        raise ValueError('refusing to write an FCS file with no events or no channels')

    # Interleave columns into the row-major order FCS requires
    data = array.array('f', bytes(4 * tot * par))
    for i, name in enumerate(names):
        data[i::par] = array.array('f', df[name].to_numpy(dtype='float32'))
    if sys.byteorder != 'little':
        data.byteswap()
    raw = data.tobytes()

    kv = {'$BYTEORD': '1,2,3,4', '$DATATYPE': 'F', '$MODE': 'L',
          '$PAR': str(par), '$TOT': str(tot)}
    originals = originals or df.attrs.get('originals', {}) or {}
    for i, name in enumerate(names, start=1):
        col_max = float(df[name].max())
        kv[f'$P{i}N'] = str(name)
        kv[f'$P{i}B'] = '32'
        kv[f'$P{i}E'] = '0,0' # float data, no log amplification
        kv[f'$P{i}R'] = str(max(1, int(math.ceil(col_max)) if math.isfinite(col_max) else 1))
        if name in originals:
            kv[f'$P{i}S'] = str(originals[name])
    for key, value in (keywords or {}).items():
        if value not in (None, ''):
            kv[str(key)] = str(value)
    kv.setdefault('$FIL', os.path.basename(file))
    for key in ('$BEGINANALYSIS', '$ENDANALYSIS', '$BEGINSTEXT', '$ENDSTEXT', '$NEXTDATA'):
        kv.setdefault(key, '0')
    kv['$BEGINDATA'] = kv['$ENDDATA'] = '0'

    # Choose a delimiter that appears in no keyword or value
    joined = ''.join(kv.keys()) + ''.join(kv.values())
    delim = next((d for d in FCS_DELIMS if d not in joined), None)
    if delim is None:
        raise ValueError('no usable TEXT delimiter: every candidate occurs in the data')

    # $BEGINDATA/$ENDDATA sit inside TEXT but describe offsets that depend on
    # TEXT's own length, so iterate until the digit counts stop changing.
    header_len = 58
    for _ in range(8):
        text = delim + delim.join(f'{k}{delim}{v}' for k, v in kv.items()) + delim
        text_start, text_end = header_len, header_len + len(text) - 1
        data_start, data_end = text_end + 1, text_end + len(raw)
        if kv['$BEGINDATA'] == str(data_start) and kv['$ENDDATA'] == str(data_end):
            break
        kv['$BEGINDATA'], kv['$ENDDATA'] = str(data_start), str(data_end)
    else:
        raise RuntimeError('FCS TEXT segment offsets failed to converge')

    def offset(n: int) -> str: # Offsets past 8 digits are reported as 0 and read from TEXT
        return str(n).rjust(8) if len(str(n)) <= 8 else '       0'

    header = ('FCS3.1' + ' ' * 4 + offset(text_start) + offset(text_end) +
              offset(data_start) + offset(data_end) + offset(0) + offset(0)).encode('ascii')
    if len(header) != header_len:
        raise RuntimeError(f'FCS header is {len(header)} bytes, expected {header_len}')

    if os.path.dirname(file):
        os.makedirs(os.path.dirname(file), exist_ok=True)
    with open(file, 'wb') as handle:
        handle.write(header)
        handle.write(text.encode('ascii', 'replace'))
        handle.write(raw)

def ncf_to_fcs(pt: str, dir: str=None, rename: dict | list | bool = True, show: bool=True) -> pd.DataFrame:
    ''' 
    ncf_to_fcs(): convert every sample holding data in a NovoCyte .ncf file to FCS 3.1 files

    Parameters:
    pt (str): path to NovoCyte .ncf file
    dir (str, optional): output directory (Default: None; an 'fcs' directory beside the .ncf)
    rename (dict | list | bool, optional): detector renaming map, list of 'OLD=NEW' strings,
                                           True for NCF_RENAME, or False for none (Default: True)
    show (bool, optional): print a line per converted sample (Default: True)

    Dependencies: os, sqlite3, pandas, ncf_samples(), ncf_events(), & fcs_write()
    '''
    if isinstance(rename, (list, tuple)): # Parse 'OLD=NEW' strings from the CLI
        parsed = {}
        for pair in rename:
            for item in str(pair).split(','):
                if not item.strip():
                    continue
                if '=' not in item:
                    raise ValueError(f'Bad rename entry {item!r}; expected OLD=NEW')
                old, new = item.split('=', 1)
                parsed[old.strip()] = new.strip()
        rename = parsed

    samples = ncf_samples(pt)
    if samples.empty:
        raise ValueError(f'No samples with acquired data in {os.path.basename(pt)}')

    con = sqlite3.connect(f'file:{pt}?mode=ro', uri=True)
    try:
        row = con.execute("SELECT FldValue FROM ENotes WHERE FldName='SW Version'").fetchone()
        sw_version = row[0] if row else None
        skipped = con.execute('''SELECT c.WellID, c.SampleName FROM SampleConfig c
                                 LEFT JOIN SampleData d ON d.SC_ID = c.SC_ID
                                 WHERE d.SD_ID IS NULL ORDER BY c.SC_ID''').fetchall()
    finally:
        con.close()
    if show and skipped:
        print(f'{len(skipped)} configured well(s) without data, skipped: ' +
              ', '.join(f'{w}/{s}' for w, s in skipped))

    source = os.path.basename(pt)
    stem = os.path.splitext(source)[0]
    dir = dir or os.path.join(os.path.dirname(os.path.abspath(pt)), 'fcs')

    files = []
    for _, sample in samples.iterrows():
        df = ncf_events(pt, int(sample['SD_ID']), parameters=sample['Parameters'], rename=rename)
        declared = int(sample['Events'] or 0)
        if declared and len(df) != declared:
            raise ValueError(f"SD_ID {sample['SD_ID']}: reassembled {len(df)} events "
                             f'but the file declares {declared}')
        date, btim = _fcs_timestamp(sample['AcquisitionTime'])
        file = os.path.join(dir, f"{_slug(stem)}_{_slug(sample['WellID'])}_{_slug(sample['SampleName'])}.fcs")
        fcs_write(df, file, keywords={
            '$CYT': 'NovoCyte', '$SRC': sample['SampleName'], '$WELLID': sample['WellID'],
            '$OP': sample['Operator'], '$VOL': sample['Volume'], '$DATE': date, '$BTIM': btim,
            'NCF_SOURCE': source, 'NCF_SD_ID': sample['SD_ID'], 'NCF_SW_VERSION': sw_version,
            'CONVERTER': 'edms.bio.fcs.ncf_to_fcs()'})
        files.append({'SD_ID': sample['SD_ID'], 'WellID': sample['WellID'],
                      'SampleName': sample['SampleName'], 'Events': len(df),
                      'PAR': df.shape[1], 'Renamed': len(df.attrs.get('originals', {})),
                      'File': file})
        if show:
            renamed = len(df.attrs.get('originals', {}))
            print(f'  wrote {os.path.basename(file)}  {len(df)} events x {df.shape[1]} params' +
                  (f'  ({renamed} channel(s) renamed)' if renamed else ''))
    if show:
        print(f'{len(files)} file(s) written to {dir}')
    return pd.DataFrame(files)

def _fcs_timestamp(timestamp) -> tuple:
    ''' 
    _fcs_timestamp(): split a NovoCyte acquisition timestamp into FCS $DATE and $BTIM

    Parameters:
    timestamp: acquisition timestamp (e.g., '2026-07-12 10:29:04')

    Dependencies: re
    '''
    if not timestamp:
        return None, None
    match = re.match(r'(\d{4})-(\d{2})-(\d{2})[ T](\d{2}:\d{2}:\d{2})', str(timestamp))
    if not match:
        return None, None
    year, month, day, time = match.groups()
    return f'{day}-{MONTHS[int(month) - 1]}-{year}', time

def _slug(text) -> str:
    ''' 
    _slug(): reduce a sample or well name to a filename-safe token

    Parameters:
    text: name to reduce (e.g., '2,H1')

    Dependencies: re
    '''
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(text or '').strip()).strip('_') or 'sample'

# Flow cytometry plots
def offset_summary(df: pd.DataFrame | str, x: str, cols: str, cols_order: list=None,
                   reps: str=None, ref: str | list=None,
                   facetx: str=None, facety: str=None, facetx_order: list=None, facety_order: list=None,
                   stat: str | bool='auto', stat_graph: str='swarm', stat_width: float=0.5,
                   stat_errorbar: str='sd', stat_axis: str=None, ref_axis: str=None,
                   x_axis: str=None, x_axis_scale: str='log', x_axis_dims: tuple | list=(0,0),
                   floor: float=1, panel: tuple=(6.5,7),
                   title: str='', title_size: int=12, file: str=None, show: bool=True,
                   dpi: int=0, transparent: bool=True, fig=None, **dist_kwargs) -> tuple:
    ''' 
    offset_summary(): draw a half-offset distribution of a flow channel beside its summary statistic
                      and, optionally, that statistic as a % difference from reference sample(s)

    Parameters:
    df (dataframe | str): tidy event dataframe (or file path); 1 row per event
    x (str): channel column name (i.e., 'PE-A')
    cols (str): sample column name; 1 baseline per value
    cols_order (list, optional): sample order, top baseline first (Default: None; appearance order)
    reps (str, optional): replicate column name; the statistic panel then holds 1 mark per replicate
                          rather than 1 over the pooled events (Default: None)
    ref (str | list, optional): cols value(s) the % difference panel reads every baseline against;
                                several are combined into 1 reference per panel & replicate
                                (Default: None; no % difference panel)
    facetx (str, optional): column name for facet columns (Default: None)
    facety (str, optional): column name for facet rows (Default: None)
    facetx_order (list, optional): facet column order (Default: None)
    facety_order (list, optional): facet row order (Default: None)
    stat (str | bool, optional): summary statistic per baseline ('gmean', 'mean', or 'auto', which is
                                 'gmean' on a log scale; Default: 'auto'). None drops the panel, which
                                 ref overrides since a % difference is of a statistic
    stat_graph (str, optional): cat() graph type the statistic is drawn as (Default: 'swarm', which
                                draws every replicate; see plot.dist() for the rest)
    stat_width (float, optional): width of each statistic panel relative to its distribution panel (Default: 0.5)
    stat_errorbar (str, optional): 'sd', 'sem' or None (Default: 'sd')
    stat_axis (str, optional): statistic panel axis name (Default: None; 'gMFI' for a geometric mean
                               of a log scaled channel, otherwise the statistic's name)
    ref_axis (str, optional): % difference panel axis name (Default: None; '% diff. vs <ref>')
    x_axis (str, optional): distribution panel axis name (Default: None; the channel column name)
    x_axis_scale (str, optional): 'log' or 'linear' (Default: 'log')
    x_axis_dims (tuple, optional): distribution panel (start, end) in channel units (Default: (0,0), which autoscales)
    floor (float, optional): clip values below this before a log scale, since a compensated channel
                             holds non-positive values that log10 cannot take (Default: 1; None skips it)
    panel (tuple, optional): (width, height) in inches of 1 distribution panel, which the statistic
                             panels & facets are sized against (Default: (6.5,7)); ignored when fig is given
    title (str, optional): figure title (Default: '')
    title_size (int, optional): title font size (Default: 12)
    file (str, optional): output file path (Default: None)
    show (bool, optional): show plot (Default: True)
    dpi (int, optional): figure dpi (Default: 0, which is plot.save_fig()'s default)
    transparent (bool, optional): save with a transparent background (Default: True)
    fig (Figure | SubFigure, optional): draw into this figure instead of making one, which ignores
                                        panel & leaves the layout, saving & showing to the caller
                                        (Default: None)
    **dist_kwargs: passed through to plot.dist() (i.e., palette_or_cmap, offset_bw_adjust, subplot_titles)

    Returns: (fig, axes) where axes holds the distribution panels; the statistic panels sit beside
             them in fig.axes, 1 distribution panel then its statistic panel(s) per facet

    Dependencies: numpy, io, plot

    Note: the figure is sized so that each distribution panel keeps `panel` whatever the panel count
    is, rather than every panel sharing a fixed figure, since a statistic panel added to a facet
    would otherwise squeeze the distribution it belongs to.
    '''
    df = io.get(pt=df) if isinstance(df, str) else df.copy()

    if cols not in df.columns:
        raise ValueError(f"cols '{cols}' is not a column! {list(df.columns)}")
    if ref is not None and not isinstance(ref, (list, tuple)):
        ref = [ref]

    # A compensated channel holds non-positive values that a log scale cannot take, & dropping
    # those events would quietly thin the very baselines that sit lowest.
    if x_axis_scale == 'log' and floor is not None:
        below = int((df[x] < floor).sum())
        if below:
            print(f'Floored {below:,} of {len(df):,} {x} events to {floor}')
        df[x] = df[x].clip(lower=floor)

    # 1 statistic panel, plus 1 more for the % difference; ref turns the statistic on by itself
    stat_panels = 0
    if ref is not None or stat not in (None, False):
        stat_panels += 1
    if ref is not None:
        stat_panels += 1
    ratios = [1] + [stat_width] * stat_panels

    # The gMFI of a log scaled channel is what a flow plot calls its geometric mean
    if stat_axis is None:
        resolved = 'gmean' if stat in (True, 'auto') and x_axis_scale == 'log' else stat
        stat_axis = 'gMFI' if resolved == 'gmean' and x_axis_scale == 'log' else ''
    if ref_axis is None and ref is not None:
        ref_axis = f"% diff. vs {ref[0] if len(ref) == 1 else 'references'}"

    # Each facet column holds a distribution panel plus its statistic panel(s), so the figure
    # grows with the panel count & the facet grid instead of the panels shrinking inside it.
    ncols = len(p._facet_values(df, facetx, order=facetx_order))
    nrows = len(p._facet_values(df, facety, order=facety_order))
    figsize = (panel[0] * sum(ratios) * ncols, panel[1] * nrows)

    fig, owns_fig = p._resolve_fig(fig, figsize)

    _, axes = p.dist(
        graph='half_offset', df=df, x=x,
        cols=cols, cols_order=cols_order,
        facetx=facetx, facety=facety, facetx_order=facetx_order, facety_order=facety_order,
        x_axis=x if x_axis is None else x_axis,
        x_axis_scale=x_axis_scale, x_axis_dims=x_axis_dims,
        offset_stat=stat, offset_stat_graph=stat_graph, offset_stat_width=stat_width,
        offset_stat_reps=reps, offset_stat_errorbar=stat_errorbar, offset_stat_axis=stat_axis,
        offset_stat_ref=ref, offset_stat_ref_axis='' if ref_axis is None else ref_axis,
        title=title, title_size=title_size,
        fig=fig, **dist_kwargs
    )

    if not owns_fig: # The caller owns the layout, saving & showing
        return fig, axes

    fig.tight_layout()
    if getattr(fig, '_suptitle', None) is not None: # Clear of the top panel, whatever its height
        top = max(ax.get_position().y1 for ax in np.ravel(axes) if ax.get_visible())
        fig._suptitle.set_y(min(0.995, top + (2.2 * title_size / 72) / fig.get_figheight()))

    p._final_save_show(fig, file=file, dpi=dpi, transparent=transparent, icon='histogram', show=show)
    return fig, axes
