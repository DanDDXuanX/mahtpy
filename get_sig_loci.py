from mahtpy import SummaryStats, MahtPlot, ColorSet
import numpy as np

import sys

add = sys.argv[1]
ofile = sys.argv[2]

ss = SummaryStats(file_input = add)

ss.get_gene(level='snp').to_csv(ofile,sep='\t',index=False)