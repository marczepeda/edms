### clone.py ###
# Author: Marc Zepeda
# Date: 2024-05-29

# Import packages
import pandas as pd
from Bio.Seq import Seq

# Supporting Methods
''' ord_form: Sigma Alrich ordering formatter
        df: Dataframe
        id: id column name
        spacer: spacer column name
        suf: suffix for oligonucleotide category
    Dependencies: pandas
'''
def ord_form(df:pd.DataFrame,id:str,seq:str,suf:str,pre:str):
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

''' tb: Designs top & bottom oligonucleotides
        df: Dataframe with sequences
        id: id column
        seq: sequence column
        t5: top oligonucleotide 5' overhang
        t3: top oligonucleotide 3' overhang
        b5: bottom oligonucleotide 5' overhang
        b3: bottom oligonucleotide 3' overhang
        tG: add 5' G to spacer if needed
        pre: prefix for ids and id column
    Dependencies: pandas, Bio.Seq
'''
def tb(df:pd.DataFrame,id:str,seq:str,t5:str,t3:str,b5:str,b3:str,tG:bool,pre:str):
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
''' sgRNAs: Design GG cloning oligonucleotides for cutting and base editing sgRNAs
        df: Dataframe with sequence information for epegRNAs
        id: id column
        spacer: spacer column (Default: Spacer_sequence)
        t5: top oligonucleotide 5' overhang
        t3: top oligonucleotide 3' overhang
        b5: bottom oligonucleotide 5' overhang (revcom)
        b3: bottom oligonucleotide 3' overhang (revcom)
        tG: add 5' G to spacer if needed (Default: True)
        order: order format
    Dependencies: pandas, top_bot(), ord_form()
'''
def sgRNAs(df:pd.DataFrame,id:str,spacer='Spacer_sequence',
          t5='CACC',t3='',b5='AAAC',b3='',
          tG=True,order=True):
    df=tb(df=df,id=id,seq=spacer,t5=t5,t3=t3,b5=b5,b3=b3,tG=tG,pre='o') # Make top and bottom oligos for spacer inserts
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='o'), # Sigma order format
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='o')]).reset_index(drop=True)
    else: return df # Original dataframe with top and bottom oligos

''' epegRNAs: Design GG cloning oligonucleotides for prime editing epegRNAs
        df: Dataframe with sequence information for epegRNAs
        id: id column
        tG: add 5' G to spacer if needed (Default: True)
        order: order format (Default: True)
        make_extension: concatenate RTT, PBS, and linker to make extension sequence (Default: True)
        }
        spacer: epegRNA spacer column (Default: Spacer_sequence)
        extension: epegRNA extension (Default: Extension_sequence)
            RTT: epegRNA reverse transcripase template column (Default: RTT_sequence)
            PBS: epegRNA primer binding site column (Default: PBS_sequence)
            linker: epegRNA linker column (Default: Linker_sequence)
        }
    Assumptions:
        epegRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
        epegRNA motif: tevoPreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA)
    Dependencies: pandas, top_bot(), ord_form()
'''
def epegRNAs(df: pd.DataFrame,id: str,tG=True, order=True,make_extension=True,
             spacer='Spacer_sequence',spacer_t5='CACC',spacer_t3='GTTTAAGAGC',spacer_b5='',spacer_b3='',
             extension='Extension_sequence',RTT='RTT_sequence',PBS='PBS_sequence',linker='Linker_sequence',extension_t5='',extension_t3='',extension_b5='CGCG',extension_b3='GCACCGACTC'):
    if make_extension==True: df[extension] = df[RTT]+df[PBS]+df[linker] # Make extension by concatenating RTT, PBS, and linker
    else: print(f'Warning: Did not make extension sequence!\nMake sure "{extension}" column includes RTT+PBS+linker for epegRNAs.')
    df=tb(df=df,id=id,seq=spacer,t5=spacer_t5,t3=spacer_t3,b5=spacer_b5,b3=spacer_b3,tG=tG,pre='ps_') # Make top and bottom oligos for spacer inserts
    df=tb(df=df,id=id,seq=extension,t5=extension_t5,t3=extension_t3,b5=extension_b5,b3=extension_b3,tG=False,pre='pe_') # Make top and bottom oligos for extension inserts
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='ps_'), # Sigma order format
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='ps_'),
                                      ord_form(df=df,id=id,seq=extension,suf='_t',pre='pe_'),
                                      ord_form(df=df,id=id,seq=extension,suf='_b',pre='pe_')]).reset_index(drop=True)
    else: return df # Original dataframe with top and bottom oligos

''' ngRNAs: Design GG cloning oligonucleotides for prime editing ngRNAs
        df: Dataframe with spacers
        id: id column
        tG: add 5' G to spacer if needed
        order: order format
        }
        spacer: ngRNA spacer column
            t5: top oligonucleotide 5' overhang
            t3: top oligonucleotide 3' overhang
            b5: bottom oligonucleotide 5' overhang
            b3: bottom oligonucleotide 3' overhang}
        }
    Assumptions:
        ngRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    Dependencies: pandas, top_bot(), ord_form()
'''
def ngRNAs(df: pd.DataFrame,id: str,tG=True, order=True,
           spacer='Spacer_sequence',ngRNA_sp_t5='CACC',ngRNA_sp_t3='GTTTAAGAGC',ngRNA_sp_b5='',ngRNA_sp_b3=''):
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
    'Homology Arm 5': 'GGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'CGGGTTCCGTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'TTCTATTCTAAGCGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'CGGGTTCCGTGGAAAGGA',
    'Reverse PCR1 Primer': 'TTCTATTCTAAGCGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 68,
    'Reverse PCR1 Primer NEB Tm': 66,
    'PCR1 Primers NEB Ta': 67,
    'Forward PCR1 Primer ID': 'MUZ141.1',
    'Reverse PCR1 Primer ID': 'MUZ141.11'
    },
    1: {
    'Subpool Number': 2,
    'Forward Barcode': 'GTTTATCGGGC',
    'Reverse Barcode': 'ACTTACTGTACC',
    'Homology Arm 5': 'GGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'GTTTATCGGGCGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'GGTACAGTAAGTGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'GTTTATCGGGCGGAAAGGA',
    'Reverse PCR1 Primer': 'GGTACAGTAAGTGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 65,
    'Reverse PCR1 Primer NEB Tm': 67,
    'PCR1 Primers NEB Ta': 66,
    'Forward PCR1 Primer ID': 'MUZ141.2',
    'Reverse PCR1 Primer ID': 'MUZ141.12'
    },
    2: {
    'Subpool Number': 3,
    'Forward Barcode': 'ACCGATGTTGAC',
    'Reverse Barcode': 'CTCGTAATAGC',
    'Homology Arm 5': 'GGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'ACCGATGTTGACGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'GCTATTACGAGGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'ACCGATGTTGACGGAAAGGA',
    'Reverse PCR1 Primer': 'GCTATTACGAGGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 66,
    'Reverse PCR1 Primer NEB Tm': 67,
    'PCR1 Primers NEB Ta': 67,
    'Forward PCR1 Primer ID': 'MUZ141.3',
    'Reverse PCR1 Primer ID': 'MUZ141.13'
    },
    3: {
    'Subpool Number': 4,
    'Forward Barcode': 'GAGGTCTTTCATGC',
    'Reverse Barcode': 'CACAACATA',
    'Homology Arm 5': 'GGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'GAGGTCTTTCATGCGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'TATGTTGTGGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'GAGGTCTTTCATGCGGAAAGGA',
    'Reverse PCR1 Primer': 'TATGTTGTGGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 67,
    'Reverse PCR1 Primer NEB Tm': 64,
    'PCR1 Primers NEB Ta': 65,
    'Forward PCR1 Primer ID': 'MUZ141.4',
    'Reverse PCR1 Primer ID': 'MUZ141.14'
    },
    4: {
    'Subpool Number': 5,
    'Forward Barcode': 'TATCCCGTGAAGCT',
    'Reverse Barcode': 'TTCGGTTAA',
    'Homology Arm 5': 'GGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'TATCCCGTGAAGCTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'TTAACCGAAGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'TATCCCGTGAAGCTGGAAAGGA',
    'Reverse PCR1 Primer': 'TTAACCGAAGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 68,
    'Reverse PCR1 Primer NEB Tm': 64,
    'PCR1 Primers NEB Ta': 65,
    'Forward PCR1 Primer ID': 'MUZ141.5',
    'Reverse PCR1 Primer ID': 'MUZ141.15'
    },
    5: {
    'Subpool Number': 6,
    'Forward Barcode': 'TAGTAGTTCAGACGC',
    'Reverse Barcode': 'ATGTACCC',
    'Homology Arm 5': 'GGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'TAGTAGTTCAGACGCGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'GGGTACATGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'TAGTAGTTCAGACGCGGAAAGGA',
    'Reverse PCR1 Primer': 'GGGTACATGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 67,
    'Reverse PCR1 Primer NEB Tm': 65,
    'PCR1 Primers NEB Ta': 66,
    'Forward PCR1 Primer ID': 'MUZ141.6',
    'Reverse PCR1 Primer ID': 'MUZ141.16'
    },
    6: {
    'Subpool Number': 7,
    'Forward Barcode': 'GGATGCATGATCTAG',
    'Reverse Barcode': 'CATCAAGC',
    'Homology Arm 5': 'GGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'GGATGCATGATCTAGGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'GCTTGATGGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'GGATGCATGATCTAGGGAAAGGA',
    'Reverse PCR1 Primer': 'GCTTGATGGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 66,
    'Reverse PCR1 Primer NEB Tm': 66,
    'PCR1 Primers NEB Ta': 67,
    'Forward PCR1 Primer ID': 'MUZ141.7',
    'Reverse PCR1 Primer ID': 'MUZ141.17'
    },
    7: {
    'Subpool Number': 8,
    'Forward Barcode': 'ATGAGGACGAATCT',
    'Reverse Barcode': 'CACCTAAAG',
    'Homology Arm 5': 'GGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'ATGAGGACGAATCTGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'CTTTAGGTGGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'ATGAGGACGAATCTGGAAAGGA',
    'Reverse PCR1 Primer': 'CTTTAGGTGGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 65,
    'Reverse PCR1 Primer NEB Tm': 65,
    'PCR1 Primers NEB Ta': 66,
    'Forward PCR1 Primer ID': 'MUZ141.8',
    'Reverse PCR1 Primer ID': 'MUZ141.18'
    },
    8: {
    'Subpool Number': 9,
    'Forward Barcode': 'GGTAGGCACG',
    'Reverse Barcode': 'TAAACTTAGAACC',
    'Homology Arm 5': 'GGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'GGTAGGCACGGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'GGTTCTAAGTTTAGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'GGTAGGCACGGGAAAGGA',
    'Reverse PCR1 Primer': 'GGTTCTAAGTTTAGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 67,
    'Reverse PCR1 Primer NEB Tm': 65,
    'PCR1 Primers NEB Ta': 66,
    'Forward PCR1 Primer ID': 'MUZ141.9',
    'Reverse PCR1 Primer ID': 'MUZ141.19'
    },
    9: {
    'Subpool Number': 10,
    'Forward Barcode': 'AGTCATGATTCAG',
    'Reverse Barcode': 'GTTGCAAGTCTAG',
    'Homology Arm 5': 'GGAAAGGACGAAACACC',
    'Homology Arm 3': 'CGCGGTTCTATCTAGTTACGCGTTAAACC',
    'Forward PCR1 Template': 'AGTCATGATTCAGGGAAAGGACGAAACACC',
    'Reverse PCR1 Template': 'CTAGACTTGCAACGGTTTAACGCGTAACTAGATAGAACCGCG',
    'Forward PCR1 Primer': 'AGTCATGATTCAGGGAAAGGA',
    'Reverse PCR1 Primer': 'CTAGACTTGCAACGGTTTAACGCGT',
    'Forward PCR1 Primer NEB Tm': 63,
    'Reverse PCR1 Primer NEB Tm': 68,
    'PCR1 Primers NEB Ta': 64,
    'Forward PCR1 Primer ID': 'MUZ141.10',
    'Reverse PCR1 Primer ID': 'MUZ141.20'
    }
}).T

''' pe_pcr2: Dataframe of PE library PCR2 primers
'''
pe_pcr2 = pd.DataFrame({
    0: {
    'Name': 'MUZ141.21',
    'Extension': 'TTTCGATTTCTTGGCTTTATATATCTT',
    'Binding': 'GTGGAAAGGACGAAACACC',
    'Description': 'PE library PCR2 FWD Primer',
    'NEB Tm': 64,
    'NEB Ta': 63
    },
    1: {
    'Name': 'MUZ141.22',
    'Extension': 'AAAAAAATTCTAGTTGGTTTAACG',
    'Binding': 'CGTAACTAGATAGAACCGCG',
    'Description': 'PE library PCR2 REV Primer',
    'NEB Tm': 62,
    'NEB Ta': 63
    }
}).T