import pandas as pd
import numpy as np

def reformat_chrom(chr_col:pd.Series)->pd.Series:
    """
    convert 'chrom' columns of SummaryStats object to int type.\n
    chrX and chrY is numbered after autosome,
    invaild value is converted to -1.

    Parameters:
    ----------
        chr_col : Series
            Series or array like, the data to convert.
    
    Returns:
    ----------
        Series
            converted chr columns
    """
    def reformat(x):
        try:
            # int type
            if type(x) is int:
                return x
            # string type
            elif type(x) is str:
                if x[0:3] == 'chr':
                    x = x[3:]
                if x == "X":
                    return 23
                if x == "Y":
                    return 24
                if x == 'M':
                    return 25
                else:
                    return int(x)
            # not supported type
            else:
                return -1
        except ValueError:
            return -1
    return np.frompyfunc(reformat,1,1)(chr_col).astype(int)

known_gene:pd.DataFrame = pd.read_csv('KnownCanonicalGene.Drop.hg38.txt', sep='\t')
known_gene['#hg38.knownCanonical.chrom'] = reformat_chrom(known_gene['#hg38.knownCanonical.chrom'])

# get mapped gene from 
def mapped(chrom,pos): 
    # if chrom == 23:
    #     chrom = 'chrX'
    # else:
    #     chrom = 'chr'+str(chrom)
    Chr_b= known_gene['#hg38.knownCanonical.chrom']==chrom
    this_chrom =  known_gene[Chr_b]
    this_chrom = this_chrom[this_chrom['hg38.kgXref.geneSymbol']!="Y_RNA"]
    Bg_b = this_chrom['hg38.knownCanonical.chromStart']<=pos
    Ed_b = this_chrom['hg38.knownCanonical.chromEnd']>=pos
    try:
        return this_chrom[Bg_b&Ed_b]['hg38.kgXref.geneSymbol'].values[0]
    except:
        return np.nan
uf_mapped = np.frompyfunc(mapped, 2, 1)

def closed(chrom,pos):
    global mapped
    mapped_gene = mapped(chrom,pos)
    if mapped_gene is np.nan:
        Chr_b= known_gene['#hg38.knownCanonical.chrom']==chrom
        this_chrom = known_gene[Chr_b]
        this_chrom = this_chrom[this_chrom['hg38.kgXref.geneSymbol']!="Y_RNA"]
        S_distance = (this_chrom['hg38.knownCanonical.chromStart'] - pos).abs()
        E_distance = (this_chrom['hg38.knownCanonical.chromEnd'] - pos).abs()
        if S_distance.min() < E_distance.min():
            K_distance:pd.Series = S_distance
            tplt = "{}>"
        else:
            K_distance:pd.Series = E_distance
            tplt = "<{}"
        closest = K_distance.idxmin()
        return tplt.format(this_chrom.loc[closest,'hg38.kgXref.geneSymbol'])
    else:
        return mapped_gene
uf_closed = np.frompyfunc(closed, 2, 1)
