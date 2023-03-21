#!/usr/bin/env python3
# coding: utf-8
# lilinxuan@genomics.cn
# summary statistics class definition

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle

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
            showgene:str            = 'gene',
            title:str               = '',
            chr_sep:int             = 20000000,
            radian:float            = 2 * np.pi,
            threshold:float         = 5e-8,
            window:int              = 500000
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
        # show gene method:
        if showgene in [False,'gene','loci']:
            if showgene == False:
                self.showgene = 'False'
            else:
                self.showgene = showgene
        # configs
        self.showgene:str       = showgene
        self.title:str          = title
        self.radian:float       = radian
        self.threshold:float    = threshold
        self.chr_sep:int        = chr_sep
        self.window:int         = window
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
            quick:bool              = True
        )->Figure:
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
        for idx, SS in enumerate(self.sumstats):
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
                self.plot_maht(axes=self.ax[idx],sumstats=plot_sumstats,quick=quick)
            # TODO 其他情况
        return self.figure
    def plot_maht(
            self,
            axes:Axes,
            sumstats:SummaryStats,
            quick:bool = True
        ):
        # hide axis
        axes.axis('off')
        # set yrange:
        ymax = max(
            -np.log10(self.threshold),
            sumstats.data['log10_pvalue'].max()
            ) * 1.2
        yunit = ymax/150
        ymin = - 5 * yunit
        # ylim = ymax - ymin
        axes.set_ylim(ymin,ymax)
        # theta
        sumstats.data['theta'] = sumstats.get_theta(chr_sep=self.chr_sep,radian=self.radian)
        # plot outframe
        # x axis:
        if sumstats.chrom > 0:
            x_min = sumstats.theta(...,sumstats.from_bp)
            x_max = sumstats.theta(...,sumstats.to_bp)
            x_sep = np.linspace(sumstats.from_bp,sumstats.to_bp,11)
            # -----
            axes.plot([x_min,x_max],[0,0],color='k')
            for x in x_sep:
                x_loc = sumstats.theta(...,x)
                axes.plot([x_loc,x_loc],[0,-2*yunit],color='k')
                axes.plot([x_loc,x_loc],[0,ymax],color='k',alpha=0.1)
                axes.text(
                    x       = x_loc,
                    y       = - 3*yunit,
                    s       = '%.2f'%(x/1e6),
                    fontsize= 14,
                    horizontalalignment = 'center',
                    verticalalignment   = 'top'
                    )
            # CHROM
            axes.text(
                x   = np.pi,
                y   = -10*yunit,
                s   = 'Chr%d(M)'%sumstats.chrom,
                fontsize = 14,
                horizontalalignment = 'center',
                verticalalignment   = 'center'
                )
        else:
            n_chrom:int = sumstats.data['chrom'].max()
            self.chrom_color = {}
            for chrom in range(1,n_chrom+1): #bar
                x_begin:float   = sumstats.theta(chrom, 0)
                x_end:float     = sumstats.theta(chrom, sumstats.chr_len[chrom])
                x_width:float   = x_end - x_begin
                x_loc:float     = (x_begin + x_end)/2
                chrom_label_rotate=self.rotate_text(x_loc)
                # color of this chromosome
                try:
                    color_this = self.chrom_color[chrom]
                except KeyError:
                    color_this = self.colorset.next()
                    self.chrom_color[chrom] = color_this
                # chromosome range
                axes.add_patch(
                    Rectangle(
                        xy      = (x_begin, -3 * yunit),
                        width   = x_width,
                        height  = 3 * yunit,
                        ec      = color_this[:3]/2,
                        color   = color_this,
                        zorder  = 1
                        )
                    )
                axes.add_patch(
                    Rectangle(
                        xy      = (x_begin, 0),
                        height  = ymax,
                        width   = x_width,
                        color   = color_this,
                        alpha   = 0.1,
                        zorder  = 1,
                    )
                    )
                # chromosome label
                chrom_label = (lambda C:str(C) if C <= sumstats.heterosome['autosome'] else sumstats.heterosome[C])(chrom)
                vert_horz = chrom<10 or chrom>sumstats.heterosome['autosome']
                axes.text(
                    x       = x_loc,
                    y       = - 4 * yunit,
                    s       = chrom_label,
                    fontsize= 14,
                    horizontalalignment = 'center',
                    verticalalignment   = 'top'
                    )
                axes.text(
                    x       = x_loc,
                    y       = 146 * yunit,
                    s       = chrom_label,
                    fontdict={
                        'size'  : {False:16,True:24}[vert_horz],
                        'weight': 'bold',
                        'color' : color_this
                        },
                    alpha   = 0.2,
                    rotation= {False:90,True:0}[vert_horz],
                    horizontalalignment = 'center',
                    verticalalignment   = 'top'
                    )
            # CHROM
            axes.text(
                x   = np.pi,
                y   = -10*yunit,
                s   = 'CHROM',
                fontsize = 14,
                horizontalalignment = 'center',
                verticalalignment   = 'center'
                )
        # y axis:
        for g in range(1,100):
            if int(ymax/g)<=10:
                break
            else:
                pass
        # =|_
        for i in range(0,int(ymax+1),g):
            axes.plot([-self.radian/120,-self.radian/300],[i,i],color='k')
            axes.text(
                x   = -self.radian/100,
                y   = i,
                s   = i,
                fontsize    = 14,
                horizontalalignment = 'right',
                verticalalignment   = 'center'
                )
        # |_
        axes.plot([-0.003*self.radian,-0.003*self.radian],[0,ymax],color='k')
        # _|
        axes.plot([self.radian*1.003,self.radian*1.003],[0,ymax],color='k')
        # _|
        axes.plot([self.radian*0.9985,self.radian*1.003],[0,0],color='k')
        # --
        axes.plot([-0.003*self.radian,1.003*self.radian],[ymax,ymax],color='k')
        axes.text(
            x   = (-0.022-0.01*np.ceil(np.log10(i))) * self.radian,
            y   = ymax/2,
            s   = '-log10(P)',
            fontsize    = 14,
            rotation    = 90,
            horizontalalignment = 'center',
            verticalalignment   = 'center'
            )
        axes.text(
            x   = self.radian/2,
            y   = ymax,
            s   = self.title+': %d'%sumstats.data['size'].median(),
            fontsize    =16,
            horizontalalignment = 'center',
            verticalalignment   = 'bottom'
            )
        # plot scatter
        if self.marker == 'scatter':
            if sumstats.chrom > 0:
                pass
            else:
                not_significant = sumstats.significant(
                    level       = 'negative',
                    threshold   = self.threshold,
                    mulitest    = self.multitest,
                    window      = self.window
                    )
                if quick:
                    not_significant = not_significant.sample(
                        n       = 50000,
                        weights = 1/not_significant['pvalue'],
                        replace = False
                        )
                genome_significant = sumstats.significant(
                    level       = 'genome',
                    threshold   = self.threshold,
                    mulitest    = self.multitest,
                    window      = self.window
                    )
                study_significant = sumstats.significant(
                    level       = 'study',
                    threshold   = self.threshold,
                    mulitest    = self.multitest,
                    window      = self.window
                    )
                for chrom in range(1,n_chrom+1):
                    this_chrom = not_significant.query("chrom==@chrom")
                    # alpha = np.frompyfunc(
                    #     lambda lgp:1-lgp/10 if lgp <=5 else 0.5,
                    #     1,1
                    # )(
                    #     this_chrom['log10_pvalue']
                    # )
                    axes.scatter(
                        x       = this_chrom['theta'],
                        y       = this_chrom['log10_pvalue'],
                        s       = 20,
                        alpha   = 1,
                        color   = self.chrom_color[chrom],
                        zorder  = 2
                        )
                axes.scatter(
                    x       = genome_significant['theta'],
                    y       = genome_significant['log10_pvalue'],
                    s       = 20,
                    alpha   = 1,
                    color   = [0.3,1,0.3,1],
                    zorder  = 2
                    )
                axes.scatter(
                    x       = study_significant['theta'],
                    y       = study_significant['log10_pvalue'],
                    s       = 30,
                    alpha   = 1,
                    color   = [1,0.3,0.3,1],
                    zorder  = 2
                    )
        # gene text
        if self.showgene != 'False':
            gene_to_show = sumstats.get_gene(
                level       = self.showgene,
                threshold   = self.threshold,
                window      = self.window
                )
            gene_to_show['theta'] = sumstats.theta(
                gene_to_show['chrom'],
                gene_to_show['chr']
            )
            # TODO: 显示基因，显示显著性阈值线

    # get the rotation angle of chromosome text label
    def rotate_text(self,theta):
        if theta <= np.pi/self.radian :
            return theta
        else:
            return theta - np.pi/self.radian