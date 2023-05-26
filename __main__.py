#!/usr/bin/env python3
# coding: utf-8
# lilinxuan@genomics.cn

# main

from core import SummaryStats, MahtPlot, ColorSet
import argparse
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Draw Manhattan plot py matplotlib!")
    # ----required----
    parser.add_argument("-i", "--input",nargs='+', help='The input summary statistics.',type=str,required=True)
    parser.add_argument("-o", "--output", help='The output plot path, format is .png',type=str,required=True)
    # ----optional----
    parser.add_argument("-t", "--threshold", help="Genome-wide significance threshold.", default=5e-8, type=np.float64)
    parser.add_argument("-c", "--cmap", help="The color schemes of mahtplot", default='default', type=str)
    parser.add_argument("-s", "--name", help="Specify the name of traits on title.", nargs='+', default='', type=str)
    parser.add_argument("--gap", help="The distance between chromosomes in plot, unit is bp.", default=20000000, type=np.int64)
    parser.add_argument("--xProtDist", help="X-axis protection distance multiplier", default=1.0, type=np.float64)
    # ----cmap----
    parser.add_argument("-n", "--ncolors", help="The number of colors in cmap.", default=2, type=np.int64)
    parser.add_argument("-p", "--posinit", help="The color pos of chr1 in plot.", default=0, type=np.int64)
    parser.add_argument("-d", "--colordist", help="The distance between adjacent colors.", default=1, type=np.int64)
    # ----region specify----
    parser.add_argument("--locus", help="Only plot region around the specified significant locus ID.", default=-1, type=np.int64)
    parser.add_argument("--chr", help="plot the specified chromosome.", default=-1, type=np.int64)
    parser.add_argument("--from-bp", help="plot region begin at this.", default=-1, type=np.int64)
    parser.add_argument("--to-bp", help="plot region end at this.", default=-1, type=np.int64)
    # ----col specify----
    parser.add_argument("--seperator", help="Specify the seperator of sumstats table.",default='\t',type=str)
    parser.add_argument("--sizecol", help="Specify a columns of input file as samplesize of GWAS summary.",default='OBS_CT',type=str)
    parser.add_argument("--chrcol", help="Specify a columns of input file as chromosome index of GWAS summary.",default='#CHROM',type=str)
    parser.add_argument("--poscol", help="Specify a columns of input file as base position of GWAS summary.",default='POS',type=str)
    parser.add_argument("--pvalcol", help="Specify a columns of input file as p-values of GWAS summary.",default='P',type=str)
    # ----flag----
    parser.add_argument("--hidegene",help="Do not show gene names.",action='store_true')
    # parse args
    args = parser.parse_args()
    # load cmap:
    if args.cmap == 'default':
        colormap = ColorSet()
    else:
        colormap = ColorSet(
            colormap    = args.cmap,
            ncolor      = args.ncolors,
            init        = args.posinit,
            dist        = args.colordist
        )
    # load sumstats
    list_of_ss = []
    for idx,file in enumerate(args.input):
        file:str = file
        try:
            trait_name = args.name[idx]
        except IndexError:
            trait_name = file.replace('\\','/').split('/')[-1]
        list_of_ss.append(
            SummaryStats(
                file_input  = file,
                seperator   = args.seperator,
                col_specify = {
                    'chrom'   : args.chrcol,
                    'pos'     : args.poscol,
                    'pvalue'  : args.pvalcol,
                    'size'    : args.sizecol,
                },
                ref_genome  = 'hg38',
                name        = trait_name
            )
        )
    # plot Manhattan plot
    if args.hidegene:
        showgene = False
    else:
        showgene = 'loci-gene'
    mhtplot = MahtPlot(
        sumstats    = list_of_ss,
        window      = 500000,
        colorset    = colormap,
        showgene    = showgene,
        geneXdist   = args.xProtDist,
        chr_sep     = args.gap
        )
    figheight = len(list_of_ss) * 6 + 2
    mhtplot.draw(
        figsize=(16,figheight),
        locus=args.locus,
        chrom=args.chr,
        from_bp=args.from_bp,
        to_bp=args.to_bp
        )
    # save fig
    mhtplot.save(args.output)