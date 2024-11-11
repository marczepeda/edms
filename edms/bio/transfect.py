### transfect.py ###
# Author: Marc Zepeda
# Date: 2024-11-07

# Import packages
import pandas as pd 
from ..gen import io as io

# PE Methods
def pe():
    epegRNAs = io.get('/Users/marczepeda/Documents/Liau_Lab/Projects/1.Miscellaneous/0.Planning/MUZ151/LC18/epegRNAs/epegRNAs_comb.csv')
    ngRNAs = io.get('/Users/marczepeda/Documents/Liau_Lab/Projects/1.Miscellaneous/0.Planning/MUZ151/LC18/ngRNAs/ngRNAs_comb.csv')
    transfection = pd.DataFrame()
    count = 0
    for (pegRNA_number,epegRNA_name) in zip(epegRNAs['pegRNA_number'],epegRNAs['Name']):
        for ngRNA_name in ngRNAs[ngRNAs['pegRNA_number']==pegRNA_number]['Name']:
            count += 1
            transfection = pd.concat([transfection,
                                    pd.DataFrame({'Transfection':[count],
                                                    'epegRNA': [epegRNA_name],
                                                    'ngRNA': [ngRNA_name]})]).reset_index(drop=True)
    

    plasmids = io.get('/Users/marczepeda/Documents/Liau_Lab/Projects/1.Miscellaneous/0.Planning/MUZ151/in/epeg_ngRNA Plasmids.csv')
    temp = pd.DataFrame()
    for (epegRNA,ngRNA) in zip(transfection['epegRNA'],transfection['ngRNA']):
        epegRNA = plasmids[plasmids['Plasmid']==epegRNA]
        ngRNA = plasmids[plasmids['Plasmid']==ngRNA]
        epegRNA = epegRNA.rename(columns={'Plasmid':'epegRNA','Colony':'epegRNA Colony','ng/uL':'epegRNA ng/uL'})
        ngRNA = ngRNA.rename(columns={'Plasmid':'ngRNA','Colony':'ngRNA Colony','ng/uL':'ngRNA ng/uL'})
        temp = pd.concat([temp,
                        pd.merge(left=epegRNA,
                                right=ngRNA,
                                on='Description')])
    transfection = pd.merge(left=transfection,
                            right=temp,
                            on=['epegRNA','ngRNA'])