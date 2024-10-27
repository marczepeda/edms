### clone.py ###
# Author: Marc Zepeda
# Date: 2024-05-29

# Import packages
import pandas as pd
from Bio.Seq import Seq

# Supporting Methods
def ord_form(df:pd.DataFrame,id:str,seq:str,suf:str,pre:str):
    ''' 
    ord_form(): Sigma Alrich ordering formatter
    
    Parameters:
    df (dataframe): pandas dataframe
    id (str): id column name
    seq (str): oligonucleotide sequence
    suf (str): suffix for oligonucleotide category
    pre (str): prefix for oligonucleotide category
    
    Dependencies: pandas
    '''
    ord = df[[(pre+id+suf),(pre+seq+suf)]]
    ord = ord.rename(columns={(pre+id+suf):'Oligo Name',(pre+seq+suf):'Sequence'})
    scale = []
    bp = []
    for s in ord['Sequence']:
        if len(s)<60: scale.append(0.025)
        else: scale.append(0.05)
        bp.append(len(s))
    ord['Scale (µmol)']=scale
    ord['bp']=bp
    return ord

def tb(df:pd.DataFrame,id:str,seq:str,t5:str,t3:str,b5:str,b3:str,tG:bool,pre:str):
    ''' 
    tb(): designs top & bottom oligonucleotides
    
    Parameters:
    df (datframe): Dataframe with sequences
    id (str): id column name
    seq (str): sequence column name
    t5 (str): top oligonucleotide 5' overhang
    t3 (str): top oligonucleotide 3' overhang
    b5 (str): bottom oligonucleotide 5' overhang
    b3 (str): bottom oligonucleotide 3' overhang
    tG (bool): add 5' G to spacer if needed
    pre (str): prefix for ids and id column

    Dependencies: pandas & Bio.Seq
    '''
    top_ids=[]
    bot_ids=[]
    top_seqs=[]
    bot_seqs=[]
    for i,s in enumerate(df[seq]):
        top_ids.append(pre+str(df.iloc[i][id])+'_t')
        bot_ids.append(pre+str(df.iloc[i][id])+'_b')
        if (tG==True)&(s[0]!='G'): top_seqs.append(t5+'G'+s+t3)
        else: top_seqs.append(t5+s+t3)
        bot_seqs.append(b5+str(Seq(s).reverse_complement()+b3))
    df[pre+id+'_t']=top_ids
    df[pre+id+'_b']=bot_ids
    df[pre+seq+'_t']=top_seqs
    df[pre+seq+'_b']=bot_seqs
    
    return df

# GG cloning
def sgRNAs(df:pd.DataFrame,id:str,spacer='Spacer_sequence',t5='CACC',t3='',b5='AAAC',b3='',tG=True,order=True):
    ''' 
    sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
    
    Parameters:
    df (dataframe): Dataframe with sequence information for epegRNAs
    id (str): id column name
    spacer (str): spacer column name (Default: Spacer_sequence)
    t5 (str): top oligonucleotide 5' overhang
    t3 (str): top oligonucleotide 3' overhang
    b5 (str): bottom oligonucleotide 5' overhang (revcom)
    b3 (str): bottom oligonucleotide 3' overhang (revcom)
    tG (bool): add 5' G to spacer if needed (Default: True)
    order (bool): order format
    
    Dependencies: pandas, top_bot(), & ord_form()
    '''
    df=tb(df=df,id=id,seq=spacer,t5=t5,t3=t3,b5=b5,b3=b3,tG=tG,pre='o') # Make top and bottom oligos for spacer inserts
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='o'), # Sigma order format
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='o')]).reset_index(drop=True)
    else: return df # Original dataframe with top and bottom oligos

def epegRNAs(df: pd.DataFrame,id: str,tG=True, order=True,make_extension=True,
             spacer='Spacer_sequence',spacer_t5='CACC',spacer_t3='GTTTAAGAGC',spacer_b5='',spacer_b3='',
             extension='Extension_sequence',RTT='RTT_sequence',PBS='PBS_sequence',linker='Linker_sequence',extension_t5='',extension_t3='',extension_b5='CGCG',extension_b3='GCACCGACTC'):
    ''' 
    epegRNAs(): design GG cloning oligonucleotides for prime editing epegRNAs
    
    Parameters:
    df (dataframe): Dataframe with sequence information for epegRNAs
    id (str): id column name
    tG (bool, optional): add 5' G to spacer if needed (Default: True)
    order (bool, optional): order format (Default: True)
    make_extension (bool, optional): concatenate RTT, PBS, and linker to make extension sequence (Default: True)
    spacer (str, optional): epegRNA spacer column name (Default: Spacer_sequence)
        _t5 (str, optional): top oligonucleotide 5' overhang
        _t3 (str, optional): top oligonucleotide 3' overhang
        _b5 (str, optional): bottom oligonucleotide 5' overhang
        _b3 (str, optional): bottom oligonucleotide 3' overhang
    extension (str, optional): epegRNA extension name (Default: Extension_sequence)
        _t5 (str, optional): top oligonucleotide 5' overhang
        _t3 (str, optional): top oligonucleotide 3' overhang
        _b5 (str, optional): bottom oligonucleotide 5' overhang
        _b3 (str, optional): bottom oligonucleotide 3' overhang
    RTT (str, optional): epegRNA reverse transcripase template column name (Default: RTT_sequence)
    PBS (str, optional): epegRNA primer binding site column name (Default: PBS_sequence)
    linker (str, optional): epegRNA linker column name(Default: Linker_sequence)
    
    Assumptions:
    1. epegRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    2. epegRNA motif: tevoPreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA)
    
    Dependencies: pandas, top_bot(), & ord_form()
    ''' 
    if make_extension==True: df[extension] = df[RTT]+df[PBS]+df[linker] # Make extension by concatenating RTT, PBS, and linker
    else: print(f'Warning: Did not make extension sequence!\nMake sure "{extension}" column includes RTT+PBS+linker for epegRNAs.')
    df=tb(df=df,id=id,seq=spacer,t5=spacer_t5,t3=spacer_t3,b5=spacer_b5,b3=spacer_b3,tG=tG,pre='ps_') # Make top and bottom oligos for spacer inserts
    df=tb(df=df,id=id,seq=extension,t5=extension_t5,t3=extension_t3,b5=extension_b5,b3=extension_b3,tG=False,pre='pe_') # Make top and bottom oligos for extension inserts
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='ps_'), # Sigma order format
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='ps_'),
                                      ord_form(df=df,id=id,seq=extension,suf='_t',pre='pe_'),
                                      ord_form(df=df,id=id,seq=extension,suf='_b',pre='pe_')]).reset_index(drop=True)
    else: return df # Original dataframe with top and bottom oligos

def ngRNAs(df: pd.DataFrame,id: str,tG=True, order=True,
           spacer='Spacer_sequence',ngRNA_sp_t5='CACC',ngRNA_sp_t3='GTTTAAGAGC',ngRNA_sp_b5='',ngRNA_sp_b3=''):
    ''' 
    ngRNAs: design GG cloning oligonucleotides for prime editing ngRNAs
    
    Parameters:
    df (dataframe): Dataframe with spacers
    id (str): id column
    tG (bool, optional): add 5' G to spacer if needed (Default: True)
    order (bool, optional): order format (Default: True)
    spacer (str, optional): ngRNA spacer column name
    ngRNA_sp_t5 (str, optional): top oligonucleotide 5' overhang
    ngRNA_sp_t3 (str, optional): top oligonucleotide 3' overhang
    ngRNA_sp_b5 (str, optional): bottom oligonucleotide 5' overhang
    ngRNA_sp_b3 (str, optional): bottom oligonucleotide 3' overhang}
    
    Assumptions:
    1. ngRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    
    Dependencies: pandas, top_bot(), & ord_form()
'''
    df=tb(df=df,id=id,seq=spacer,t5=ngRNA_sp_t5,t3=ngRNA_sp_t3,b5=ngRNA_sp_b5,b3=ngRNA_sp_b3,tG=tG,pre='ns_') # Make top and bottom oligos for spacer inserts
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='ns_'), # Sigma order format
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='ns_')]).reset_index(drop=True)
    
    else: return df # Original dataframe with top and bottom oligos

# Gibson cloning
''' pe_pcr1: Dataframe of PE library PCR1 primers
'''
pe_pcr1 = pd.DataFrame({
    0: {
    'Subpool Number': 1,
    'Forward Barcode': 'CGGGTTCCGT',
    'Reverse Barcode': 'GCTTAGAATAGAA',
    'Homology Arm 5': 'GTGGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'CGGGTTCCGTGTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'TTCTATTCTAAGCGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'CGGGTTCCGTGTGGAAAG',
    'Reverse PCR1 Primer': 'TTCTATTCTAAGCGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 66,
    'Reverse PCR1 Primer NEB Tm': 66,
    'PCR1 Primers NEB Ta': 67,
    'Forward PCR1 Primer IDT Tm': 62.5,
    'Reverse PCR1 Primer IDT Tm': 64.7,
    'Forward PCR1 Primer IDT Th': 51.4,
    'Reverse PCR1 Primer IDT Th': 51.4,
    'Forward PCR1 Primer ID': 'fpMUZ141.1',
    'Reverse PCR1 Primer ID': 'rpMUZ141.1'
    },
    1: {
    'Subpool Number': 2,
    'Forward Barcode': 'GTTTATCGGGC',
    'Reverse Barcode': 'ACTTACTGTACC',
    'Homology Arm 5': 'GTGGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'GTTTATCGGGCGTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'GGTACAGTAAGTGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'GTTTATCGGGCGTGGAAAG',
    'Reverse PCR1 Primer': 'GGTACAGTAAGTGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 64,
    'Reverse PCR1 Primer NEB Tm': 67,
    'PCR1 Primers NEB Ta': 65,
    'Forward PCR1 Primer IDT Tm': 60.8,
    'Reverse PCR1 Primer IDT Tm': 64.7,
    'Forward PCR1 Primer IDT Th': -4.3,
    'Reverse PCR1 Primer IDT Th': -4.3,
    'Forward PCR1 Primer ID': 'fpMUZ141.2',
    'Reverse PCR1 Primer ID': 'rpMUZ141.2'
    },
    2: {
    'Subpool Number': 3,
    'Forward Barcode': 'ACCGATGTTGAC',
    'Reverse Barcode': 'CTCGTAATAGC',
    'Homology Arm 5': 'GTGGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'ACCGATGTTGACGTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'GCTATTACGAGGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'ACCGATGTTGACGTGGAAAG',
    'Reverse PCR1 Primer': 'GCTATTACGAGGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 65,
    'Reverse PCR1 Primer NEB Tm': 67,
    'PCR1 Primers NEB Ta': 66,
    'Forward PCR1 Primer IDT Tm': 62.5,
    'Reverse PCR1 Primer IDT Tm': 64.7,
    'Forward PCR1 Primer IDT Th': 36.8,
    'Reverse PCR1 Primer IDT Th': 36.8,
    'Forward PCR1 Primer ID': 'fpMUZ141.3',
    'Reverse PCR1 Primer ID': 'rpMUZ141.3'
    },
    3: {
    'Subpool Number': 4,
    'Forward Barcode': 'GAGGTCTTTCATGC',
    'Reverse Barcode': 'CACAACATA',
    'Homology Arm 5': 'GTGGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'GAGGTCTTTCATGCGTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'TATGTTGTGGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'GAGGTCTTTCATGCGTGGAAAG',
    'Reverse PCR1 Primer': 'TATGTTGTGGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 66,
    'Reverse PCR1 Primer NEB Tm': 64,
    'PCR1 Primers NEB Ta': 65,
    'Forward PCR1 Primer IDT Tm': 63.4,
    'Reverse PCR1 Primer IDT Tm': 62.7,
    'Forward PCR1 Primer IDT Th': 45.3,
    'Reverse PCR1 Primer IDT Th': 45.3,
    'Forward PCR1 Primer ID': 'fpMUZ141.4',
    'Reverse PCR1 Primer ID': 'rpMUZ141.4'
    },
    4: {
    'Subpool Number': 5,
    'Forward Barcode': 'TAGTAGTTCAGACGC',
    'Reverse Barcode': 'ATGTACCC',
    'Homology Arm 5': 'GTGGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'TAGTAGTTCAGACGCGTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'GGGTACATGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'TAGTAGTTCAGACGCGTGGAAAG',
    'Reverse PCR1 Primer': 'GGGTACATGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 66,
    'Reverse PCR1 Primer NEB Tm': 65,
    'PCR1 Primers NEB Ta': 66,
    'Forward PCR1 Primer IDT Tm': 63.8,
    'Reverse PCR1 Primer IDT Tm': 62.6,
    'Forward PCR1 Primer IDT Th': 18.5,
    'Reverse PCR1 Primer IDT Th': 18.5,
    'Forward PCR1 Primer ID': 'fpMUZ141.5',
    'Reverse PCR1 Primer ID': 'rpMUZ141.5'
    },
    5: {
    'Subpool Number': 6,
    'Forward Barcode': 'GGATGCATGATCTAG',
    'Reverse Barcode': 'CATCAAGC',
    'Homology Arm 5': 'GTGGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'GGATGCATGATCTAGGTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'GCTTGATGGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'GGATGCATGATCTAGGTGGAAAG',
    'Reverse PCR1 Primer': 'GCTTGATGGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 65,
    'Reverse PCR1 Primer NEB Tm': 66,
    'PCR1 Primers NEB Ta': 66,
    'Forward PCR1 Primer IDT Tm': 62.7,
    'Reverse PCR1 Primer IDT Tm': 63.4,
    'Forward PCR1 Primer IDT Th': 35.3,
    'Reverse PCR1 Primer IDT Th': 35.3,
    'Forward PCR1 Primer ID': 'fpMUZ141.6',
    'Reverse PCR1 Primer ID': 'rpMUZ141.6'
    },
    6: {
    'Subpool Number': 7,
    'Forward Barcode': 'ATGAGGACGAATCT',
    'Reverse Barcode': 'CACCTAAAG',
    'Homology Arm 5': 'GTGGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'ATGAGGACGAATCTGTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'CTTTAGGTGGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'ATGAGGACGAATCTGTGGAAAG',
    'Reverse PCR1 Primer': 'CTTTAGGTGGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 64,
    'Reverse PCR1 Primer NEB Tm': 65,
    'PCR1 Primers NEB Ta': 65,
    'Forward PCR1 Primer IDT Tm': 62.2,
    'Reverse PCR1 Primer IDT Tm': 62.8,
    'Forward PCR1 Primer IDT Th': 24.0,
    'Reverse PCR1 Primer IDT Th': 24.0,
    'Forward PCR1 Primer ID': 'fpMUZ141.7',
    'Reverse PCR1 Primer ID': 'rpMUZ141.7'
    },
    7: {
    'Subpool Number': 8,
    'Forward Barcode': 'GGTAGGCACG',
    'Reverse Barcode': 'TAAACTTAGAACC',
    'Homology Arm 5': 'GTGGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'GGTAGGCACGGTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'GGTTCTAAGTTTAGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'GGTAGGCACGGTGGAAAG',
    'Reverse PCR1 Primer': 'GGTTCTAAGTTTAGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 66,
    'Reverse PCR1 Primer NEB Tm': 65,
    'PCR1 Primers NEB Ta': 66,
    'Forward PCR1 Primer IDT Tm': 62.0,
    'Reverse PCR1 Primer IDT Tm': 63.4,
    'Forward PCR1 Primer IDT Th': 13.2,
    'Reverse PCR1 Primer IDT Th': 13.2,
    'Forward PCR1 Primer ID': 'fpMUZ141.8',
    'Reverse PCR1 Primer ID': 'rpMUZ141.8'
    },
    8: {
    'Subpool Number': 9,
    'Forward Barcode': 'AGTCATGATTCAG',
    'Reverse Barcode': 'GTTGCAAGTCTAG',
    'Homology Arm 5': 'GTGGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'AGTCATGATTCAGGTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'CTAGACTTGCAACGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'AGTCATGATTCAGGTGGAAAG',
    'Reverse PCR1 Primer': 'CTAGACTTGCAACGGTTTAACGC',
    'Forward PCR1 Primer NEB Tm': 62,
    'Reverse PCR1 Primer NEB Tm': 66,
    'PCR1 Primers NEB Ta': 63,
    'Forward PCR1 Primer IDT Tm': 60.1,
    'Reverse PCR1 Primer IDT Tm': 63.9,
    'Forward PCR1 Primer IDT Th': 36.5,
    'Reverse PCR1 Primer IDT Th': 36.5,
    'Forward PCR1 Primer ID': 'fpMUZ141.9',
    'Reverse PCR1 Primer ID': 'rpMUZ141.9'
    }
}).T

''' pe_pcr2: Dataframe of PE library PCR2 primers
'''
pe_pcr2 = pd.DataFrame({
    0: {
    'Name': 'fpMUZ141.0',
    'Extension': 'GCTTTATATATCTT',
    'Binding': 'GTGGAAAGGACGAAACACC',
    'Description': 'PE library PCR2 FWD Primer',
    'NEB Tm': 64,
    'NEB Ta': 63
    },
    1: {
    'Name': 'rpMUZ141.0',
    'Extension': '',
    'Binding': 'GGTTTAACGCGTAACTAGATAGAA',
    'Description': 'PE library PCR2 REV Primer',
    'NEB Tm': 62,
    'NEB Ta': 63
    }
}).T

def pe_twist_oligos(df: pd.DataFrame,tG=True, make_extension=True,spacer='Spacer_sequence',scaffold='Scaffold_sequence',
                    extension='Extension_sequence',RTT='RTT_sequence',PBS='PBS_sequence',linker='Linker_sequence'):
    ''' 
    pe_twist_oligos: makes prime editing twist olionucleotides
    
    Parameters:
    df (dataframe): Dataframe with sequence information for epegRNAs
    tG (bool, optional): add 5' G to spacer if needed (Default: True)
    make_extension (bool, optional): concatenate RTT, PBS, and linker to make extension sequence (Default: True)
    spacer (str, optional): epegRNA spacer column name (Default: Spacer_sequence)
    scaffold (str, optional): epegRNA scaffold column name (Default: Scaffold_sequence)
    extension (str, optional): epegRNA extension name (Default: Extension_sequence)
    RTT (str, optional): epegRNA reverse transcripase template column name (Default: RTT_sequence)
    PBS (str, optional): epegRNA primer binding site column name (Default: PBS_sequence)
    linker (str, optional): epegRNA linker column name (Default: Linker_sequence)
    
    Assumptions:
    1. epegRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    2. epegRNA motif: tevoPreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA)
    
    Dependencies: pandas, pe_pcr1, & pe_pcr2
    '''
    # Make extension by concatenating RTT, PBS, and linker
    if make_extension==True: df[extension] = df[RTT]+df[PBS]+df[linker]
    else: print(f'Warning: Did not make extension sequence!\nMake sure "{extension}" column includes RTT+PBS+linker for epegRNAs.')

    # Assign subpool barcodes to different PBS lengths and strand directions
    pbs_lengths = sorted(df['PBS_length'].value_counts().keys())
    strands = sorted(df['Strand'].value_counts().keys())
    pbs_lengths_ls = pbs_lengths*len(strands)
    strands_ls = strands*len(pbs_lengths)
    pcr1_barcodes = pe_pcr1.iloc[:len(pbs_lengths_ls)]
    pcr1_barcodes['PBS_length']=pbs_lengths_ls
    pcr1_barcodes['Strand']=strands_ls
    df = pd.merge(left=df,right=pcr1_barcodes,on=['PBS_length','Strand'])

    # Make twist oligo & determine length
    df['First_spacer_nucleotide']=[s[0] for s in df[spacer]]
    if tG==True: # Append 5'G to spacer if not already present
        df['Twist_oligo']=[df.iloc[i]['Forward Barcode']+df.iloc[i]['Homology Arm 5']+df.iloc[i][spacer]+
                           df.iloc[i][scaffold]+df.iloc[i][extension]+df.iloc[i]['Homology Arm 3']+
                           df.iloc[i]['Reverse Barcode'] if df.iloc[i]['First_spacer_nucleotide']=='G' else 
                           df.iloc[i]['Forward Barcode']+df.iloc[i]['Homology Arm 5']+'G'+df.iloc[i][spacer]+
                           df.iloc[i][scaffold]+df.iloc[i][extension]+df.iloc[i]['Homology Arm 3']+
                           df.iloc[i]['Reverse Barcode'] for i in range(len(df[spacer]))]
    else: # Do not append 5'G to spacer if not already present
        df['Twist_oligo']=[df.iloc[i]['Forward Barcode']+df.iloc[i]['Homology Arm 5']+df.iloc[i][spacer]+
                           df.iloc[i][scaffold]+df.iloc[i][extension]+df.iloc[i]['Homology Arm 3']+
                           df.iloc[i]['Reverse Barcode'] for i in range(len(df[spacer]))]
    df['twist_oligo_length']=[len(twist) for twist in df['Twist_oligo']]

    return df