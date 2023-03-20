#!/usr/bin/env python3
# coding: utf-8
# lilinxuan@genomics.cn
# summary statistics class definition

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from .SummaryStats import SummaryStats,SlicedSumStats
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
                # first summary, get from specified range
                plot_sumstats = SS.slice(
                    chrom = chrom,
                    from_bp = from_bp,
                    to_bp = to_bp,
                    locus = locus,
                    threshold = self.threshold,
                    window = self.window
                )
                # recode the return range
                chrom = plot_sumstats.chrom
                from_bp = plot_sumstats.from_bp
                to_bp = plot_sumstats.to_bp
            else:
                plot_sumstats = SS.slice(
                    chrom = chrom,
                    from_bp = from_bp,
                    to_bp = to_bp,
                    locus = -1,
                    threshold = self.threshold,
                    window = self.window
                )
            if self.layout in ['vertical','horizontal'] and self.style == 'classic':
                self.chr_color = {}
                self.plot_maht(axes=self.ax[idx],sumstats=plot_sumstats)
            # TODO 其他情况
    def plot_maht(
            self,
            axes:Axes,
            sumstats:SummaryStats,
        ):
        # hide axis
        axes.axis('off')
        # set yrange:
        ymin = - sumstats.data['log10_pavlue'].max() * 0.03
        ymax = max(
            -np.log10(self.threshold),
            sumstats.data['log10_pavlue'].max()
            ) * 1.2
        ylim = ymax - ymin
        # theta
        theta = sumstats.get_theta(chr_sep=self.chr_sep,radian=self.radian)
        # plot outframe
        # x axis:
        if len(sumstats.chrom > 0):
            xmin = sumstats.from_bp
            xmax = sumstats.to_bp
        else:
            n_chrom:int = sumstats.data['chrom'].max()
            self.chrom_color = {}
            for chrom in range(1,n_chrom+1): #bar
                x_loc:float = sumstats.theta(chrom, sumstats.chr_len[chrom]/2)
                chrom_label_rotate=self.rotate_text(x_loc)
                # color of this chromosome
                try:
                    color_this = self.chrom_color[chrom]
                except KeyError:
                    color_this = self.colorset.next()
                    self.chrom_color[chrom] = color_this
                # chromosome range
                axes.bar(
                    x       = x_loc,
                    height  = ylim * 3/150,
                    width   = sumstats.chr_len[chrom]/sumstats.chr_len_total*self.radian,
                    bottom  = -ylim/50,
                    ec      = color_this[:3]/2,
                    color   = color_this,
                    zorder  = 3
                    )
                axes.bar(
                    x       = x_loc,
                    height  = ymax,
                    width   = sumstats.chr_len[chrom]/sumstats.chr_len_total*self.radian,
                    bottom  = 0,
                    color   = color_this,
                    zorder  = 3,
                    alpha   = 0.1
                    )
                # chromosome label
                chrom_label = (lambda C:str(C) if C < sumstats.heterosome['autosome'] else sumstats.heterosome[C])(chrom)
                vert_horz = chrom<10 or chrom>sumstats.heterosome['autosome']
                axes.text(
                    x       = x_loc,
                    y       = - ylim * 4/150,
                    s       = 'chr'+chrom_label,
                    fontsize= 14,
                    horizontalalignment = 'center',
                    verticalalignment   = 'top'
                    )
                axes.text(
                    x       = x_loc,
                    y       = ymax-ylim*4/150,
                    s       = chrom_label,
                    fontdict={
                        'size'  : {False:16,True:24}[vert_horz],
                        'weight': 'bold',
                        'color' : color_this
                        },
                    alpha   = 0.2
                    rotation= {False:90,True:0}[vert_horz]
                    horizontalalignment = 'center',
                    verticalalignment   = 'top'
                )
    # get the rotation angle of chromosome text label
    def rotate_text(self,theta):
        if theta <= np.pi/self.radian :
            return theta
        else:
            return theta - np.pi/self.radian