''' 
Module: ncbi.py
Author: Marc Zepeda
Created: 2024-12-25
Description: National Center of Biotechnology Information

Usage:
[NCBI Information]
- get_info(): Get NCBI information
- set_info(): Set NCBI information
'''

# Import packages
from Bio import Entrez
import pandas as pd

from ..gen import tidy as t
from .. import config

# NCBI Information
def get_info():
    '''
    get_info(): Get NCBI information
    
    Dependencies: config
    '''
    # Get NCBI information
    info_ncbi = config.get_info(id='NCBI')

    # Load NCBI information
    if 'email' in list(info_ncbi.keys()): Entrez.email = info_ncbi['email']
    if 'tool' in list(info_ncbi.keys()): Entrez.tool = info_ncbi['tool']
    if 'api_key' in list(info_ncbi.keys()): Entrez.api_key = info_ncbi['api_key']
    if 'max_tries' in list(info_ncbi.keys()): Entrez.max_tries = info_ncbi['max_tries']
    if 'sleep_between_tries' in list(info_ncbi.keys()): Entrez.sleep_between_tries = info_ncbi['sleep_between_tries']

def set_info(email:str=None, tool:str=None, api_key:str=None, max_tries:int=None, sleep_between_tries:int=None):
    '''
    set_info(): Set NCBI information
    
    Parameters:
    email (str, optional): email address (Default: None)
    tool (str, optional): tool (Default: None, 'biopython')
    api_key (str, optional): Personal API key from NCBI (Default: None)
        If not set, only 3 queries per second are allowed.
        10 queries per seconds otherwise with a valid API key.
    max_tries (int, optional): configures how many times failed requests will be automatically retried on error (Default: None, '3')
    sleep_between_tries (int, optional): delay, in seconds, before retrying a request on error (Default: None, '15')

    Dependencies: get_ncbi()
    '''
    # Get previous NCBI information
    info_ncbi = config.get_info(id='NCBI')
    if info_ncbi is None: info_ncbi=dict()

    # Update NCBI information based on parameters
    if email: info_ncbi['email']=email
    if tool: info_ncbi['tool']=tool
    if api_key: info_ncbi['api_key']=api_key
    if max_tries: info_ncbi['max_tries']=max_tries
    if sleep_between_tries: info_ncbi['sleep_between_tries']=sleep_between_tries
    config.set_info(id='NCBI',info=info_ncbi)
    
    # Load NCBI information
    if 'email' in list(info_ncbi.keys()): Entrez.email = info_ncbi['email']
    if 'tool' in list(info_ncbi.keys()): Entrez.tool = info_ncbi['tool']
    if 'api_key' in list(info_ncbi.keys()): Entrez.api_key = info_ncbi['api_key']
    if 'max_tries' in list(info_ncbi.keys()): Entrez.max_tries = info_ncbi['max_tries']
    if 'sleep_between_tries' in list(info_ncbi.keys()): Entrez.sleep_between_tries = info_ncbi['sleep_between_tries']

def data(func: str, **kwargs):
    '''
    data(): retrieve data from NCBI

    Parameters:
    func (str): Entrez function name
    - efetch: Retrieves records in the requested format from a list of one or more primary IDs or from the user's environment
        - db (str): Database to fetch from.
	    - id (str): List of IDs to retrieve (comma-separated).
        - rettype (str): Data type to return (e.g., gb, fasta).
	    - retmode (str, 'xml'): Format of the output (e.g., text, xml).
	    - seq_start (int): Starting sequence position (for sequence records).
	    - seq_stop (int): Ending sequence position (for sequence records).
	    - strand (int): Strand of DNA (1 for plus, 2 for minus).
    - epost: Posts a file containing a list of primary IDs for future use in the user's environment to use with subsequent search strategies
        - db (str): Database where IDs belong.
	    - id (str): List of IDs to post.
	    - query_key (str): Key for the query results.
	    - WebEnv (str): Web environment (history) string.
    - esearch: Searches and retrieves primary IDs (for use in EFetch, ELink, and ESummary) and term translations and optionally retains results 
      for future use in the user's environment.
	    - db (str): Database to search.
	    - term (str): Search query.
	    - retmax (int, optional): Maximum number of IDs to return.
	    - retstart (int, optional): Starting point for results.
	    - usehistory (str, optional): Whether to store results on the NCBI History server (y or n).
	    - idtype (str, optional): Type of IDs to return (e.g., acc for accession numbers).
    - elink: Checks for the existence of an external or Related Articles link from a list of one or more primary IDs.  Retrieves primary IDs
      and relevancy scores for links to Entrez databases or Related Articles;  creates a hyperlink to the primary LinkOut provider for a specific 
      ID and database, or lists LinkOut URLs and Attributes for multiple IDs.
        - dbfrom (str): Source database.
	    - db (str): Target database.
	    - id (str, comma-seperated): List of IDs to link from.
	    - linkname (str): Specific link name for database relationships (e.g., pubmed_pubmed).
    - einfo: Provides field index term counts, last update, and available links for each database.
        - db: (str, Optional) Database name to query.
    - esummary: Retrieves document summaries from a list of primary IDs or from the user's environment.
        - db (str): Database to query.
	    - id (str, comma-seperated): List of IDs to summarize.
	    - retmode (str): Format of the summary (e.g., xml, json).
    - egquery: Provides Entrez database counts in XML for a single search using Global Query.
        - term (str): Search query.
    - espell: Retrieves spelling suggestions.
        - db (str, Required): Database to query for spelling suggestions.
        - term (str, Required): The query string or search terms for which you want spelling suggestions.
    - ecitmatch: Retrieves PubMed IDs (PMIDs) that correspond to a set of input citation strings.
        - db (str, Required): Must be "pubmed" (only works with PubMed).
	    - bdata (str, Required): Citation data in a specific format (pipe-delimited). Each citation must include at least a journal title, year, volume, and page number
    
    NCBI Databases: 'pubmed', 'protein', 'nuccore', 'ipg', 'nucleotide', 'structure', 'genome', 'annotinfo', 'assembly', 'bioproject', 'biosample', 
    'blastdbinfo', 'books', 'cdd', 'clinvar', 'gap', 'gapplus', 'grasp', 'dbvar', 'gene', 'gds', 'geoprofiles', 'medgen', 'mesh', 'nlmcatalog', 
    'omim', 'orgtrack', 'pmc', 'popset', 'proteinclusters', 'pcassay', 'protfam', 'pccompound', 'pcsubstance', 'seqannot', 'snp', 'sra', 'taxonomy', 
    'biocollections', and 'gtr'

    Dependencies: Entrez, tidy
    '''
    # Load NCBI credentials
    get_info()

    # Compute Function
    if func=='efetch':
        # Filter kwargs
        kwargs = t.filter_kwargs(keywords=['db','id','rettype','seq_start','seq_stop','strand'],kwargs=kwargs)
        kwargs['retmode'] = 'xml'
        
        # Interact with NCBI API
        handle = Entrez.efetch(**kwargs)
        records = Entrez.parse(handle)
        handle.close()

        # Return result
        df = pd.DataFrame()
        for record in records: df = pd.concat([df,pd.Series(record)],axis=1)
        df = df.T
        print(f"efetch:\n{df}")
        return df
        
    elif func=='epost': # Not useful...
        # Filter kwargs
        kwargs = t.filter_kwargs(keywords=['db','id','query_key','WebEnv'],kwargs=kwargs)

        # Interact with NCBI API
        handle = Entrez.epost(**kwargs)
        record = handle.read()
        handle.close()

        # Return result
        print(f"epost:\n{record}")
        return record
    
    elif func=='esearch':
        # Filter kwargs
        kwargs = t.filter_kwargs(keywords=['db','term','retmax','retstart','usehistory','idtype'],kwargs=kwargs)
        
        # Interact with NCBI API
        handle = Entrez.esearch(**kwargs)
        record = Entrez.read(handle)
        handle.close()

        # Return result
        print(f"esearch:\n{record}")
        return record
        
    elif func=='elink':
        # Filter kwargs
        kwargs = t.filter_kwargs(keywords=['dbfrom','db','id','linkname'],kwargs=kwargs)
        
        # Interact with NCBI API
        handle = Entrez.elink(**kwargs)
        record = Entrez.read(handle)
        handle.close()

        # Return result
        dc = t.comb_dcs(record[0]['LinkSetDb'][0]['Link'])
        dc['To_ID'] = dc.pop('Id')
        dc['From_ID'] = [', '.join(record[0]['IdList'])]*len(dc['To_ID'])
        dc['To_DB'] = [record[0]['LinkSetDb'][0]['DbTo']]*len(dc['To_ID'])
        dc['From_DB'] = [record[0]['DbFrom']]*len(dc['To_ID'])
        dc['Link_Name'] = [record[0]['LinkSetDb'][0]['LinkName']]*len(dc['To_ID'])
        df = pd.DataFrame(dc)
        print(f"elink:\n{df}")
        return df
    
    elif func=='einfo':
        # Filter kwargs
        kwargs = t.filter_kwargs(keywords=['db'],kwargs=kwargs)

        # Interact with NCBI API
        handle = Entrez.einfo(**kwargs)
        record = Entrez.read(handle)
        handle.close()

        # Return result
        print(f"einfo:\n{record}")
        return record
    
    elif func=='esummary':
        # Filter kwargs
        kwargs = t.filter_kwargs(keywords=['db','id','retmode'],kwargs=kwargs)
        
        # Interact with NCBI API
        handle = Entrez.esummary(**kwargs)
        record = Entrez.read(handle)
        handle.close()

        # Return result
        df = pd.DataFrame(record)
        print(f"esummary:\n{df}")
        return df

    elif func=='egquery': # Does not work
        # Filter kwargs
        kwargs = t.filter_kwargs(keywords=['term'],kwargs=kwargs)

        # Return result
        print(f"egquery: does not work unfortunately")
        return None

    elif func=='espell':
        # Filter kwargs
        kwargs = t.filter_kwargs(keywords=['db','term'],kwargs=kwargs)

        # Interact with NCBI API
        record = Entrez.read(Entrez.espell(**kwargs))

        # Return result
        df = pd.DataFrame([record])
        print(f"espell:\n{df}")
        return df

    elif func=='ecitmatch':
        # Filter kwargs
        kwargs = t.filter_kwargs(keywords=['db','bdata'],kwargs=kwargs)

        # Interact with NCBI API
        handle = Entrez.ecitmatch(**kwargs)
        record = handle.read().strip().split("|")
        handle.close()

        # Return result
        print(f"ecitmatch:\n{record}")
        return record

    else: TypeError(f"{func} is not efetch, epost, esearch, elink, einfo, esummary, egquery, espell, and ecitmatch.")