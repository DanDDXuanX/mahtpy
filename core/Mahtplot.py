#!/usr/bin/env python3
# coding: utf-8
# lilinxuan@genomics.cn
# summary statistics class definition

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from .SummaryStats import SummaryStats
from .Cmap import ColorSet

import warnings

class MahtplotError(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class MahtPlot:
    def __init__(
            self,
            sumstats:SummaryStats,
            colorset:ColorSet       = ColorSet(),
            style:str               = 'classic',
            layout:str              = 'horizontal',
            marker:str              = 'scatter',
            showgene:bool           = True,
            title:str               = '',
            chr_sep:int             = 20000000,
            radian:float            = 2 * np.pi,
            threshold:float         = 5e-8,
            ) -> None:
        # summary
        if type(sumstats) is SummaryStats:
            self.sumstats:list = [sumstats]
            self.multitest:int = 1
        elif type(sumstats) in [list,tuple,np.ndarray]:
            # check if all element of list is SummaryStats
            for ss in sumstats:
                if type(ss) is not SummaryStats:
                    raise MahtplotError("Invalid sumstats type: {}".format(type(ss)))
                if ss.ref_genome != last.ref_genome:
                    raise MahtplotError(
                        "Reference Genome ({} and {}) of input SumStats did not match."
                        .format(ss.ref_genome,last.ref_genome)
                        )
                last = ss
            self.sumstats:list = list(sumstats)
            self.multitest:int = len(sumstats)
        else:
            raise MahtplotError('Invalid sumstats argument.')
        # init the rcParams
        plt.rcParams['figure.subplot.left']=0.1
        plt.rcParams['figure.subplot.right']=0.9
        plt.rcParams['figure.subplot.top']=0.9
        plt.rcParams['figure.subplot.bottom']=0.1
        # colorset
        self.colorset:ColorSet = colorset
        # style and layout
        if style in ['classic','overlap','symmetric']:
            self.style:str = style
        else:
            raise MahtplotError("Invalid mahtplot style '{}'".format(style))
        if layout in ['vertical','horizontal','cavern','firework']:
            self.layout:str = layout
        else:
            raise MahtplotError("Invalid mahtplot layout '{}'".format(layout))
        if marker in ['scatter','lolipop']:
            self.marker:str = marker
        else:
            raise MahtplotError("Invalid mahtplot marker '{}'".format(layout))
        # configs
        self.showgene:bool      = showgene
        self.title:str          = title
        self.radian:float       = radian
        self.threshold:float    = threshold
        self.chr_sep:int        = chr_sep
        # MahtPlot class is not suitable to show overmany summary
        if self.multitest >= 5:
            warnings.warn('''
            MahtPlot class is not suitable to display more than 4 summary, 
            Use CircosMPlot instead.
            ''')
        elif self.style == 'symmetric' and self.multitest > 2:
            raise MahtplotError("Mahtplot layout 'symmetric' can only display 2 or less SumStats")

    def draw(
            self,
            figsize:tuple           = (16,8),
            chrom:int               = -1,
            from_bp:int             = -1,
            to_bp:int               = -1,
            locus:int               = -1,
        ):
        # figure and axes
        self.figure:Figure = plt.figure(figsize=figsize)
        # if only 1 SumStats, or use circular layout
        if self.multitest == 1 or self.layout in ['cavern','firework'] or self.style in ['symmetric','overlap']:
            self.ax:list = [self.figure.subplots(1,1)]
        # if mulit SumStats
        else:
            if self.layout == 'vertical':
                self.ax:list = self.figure.subplots(1, self.multitest)
            else:
                self.ax:list = self.figure.subplots(self.multitest, 1)
        # TODO: 根据 chrom 和 locus 的指定，取SumStats的subset，再考虑画图。
        for idx, SS in self.sumstats:
            SS:SummaryStats
            if idx == 0:
                plot_sumstats, chrange = SS.get_theta(
                    chr_sep = self.chr_sep,
                    radian  = self.radian,
                    chrom = chrom,
                    from_bp = from_bp,
                    to_bp = to_bp,
                    locus = locus,
                    threshold = self.threshold,
                    window = self.window
                )
    def plot_maht(
            self,
            axes:Axes,
            sumstats:SummaryStats,
        ):
        # hide axis
        axes.axis('off')
        # set yrange:
        ymin = -self.sumstats['log10_pavlue'].max()*0.03
        ymax = max(
            -np.log10(self.threshold),
            self.sumstats['log10_pavlue'].max()
            ) * 1.2