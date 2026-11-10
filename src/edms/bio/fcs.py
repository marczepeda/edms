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
import pandas as pd

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
