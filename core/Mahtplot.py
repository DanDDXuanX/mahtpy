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
            chr_sep:int             = 20000000,
            radian:float            = 2 * np.pi,
            threshold:float         = 5e-8,
            window:int              = 500000,
            geneXdist:float         = 1.0
            ) -> None:
        # summary
        if type(sumstats) is SummaryStats:
            self.sumstats:list = [sumstats]
            self.multitest:int = 1
        elif type(sumstats) in [list,tuple,np.ndarray]:
            # check if all element of list is SummaryStats
            last = sumstats[0]
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
        if showgene in [False,'gene','loci','loci-gene']:
            if showgene == False:
                self.showgene = 'False'
            else:
                self.showgene = showgene
        # configs
        self.radian:float       = radian
        self.threshold:float    = threshold
        self.chr_sep:int        = chr_sep
        self.window:int         = window
        self.geneXdist:int      = geneXdist
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
        ):
        # figsize
        self.width  = figsize[0]
        self.height = figsize[1]
        self.WHR    = self.width/self.height
        self.xzoom  = self.width/16
        self.yzoom  = self.height/8
        # init the colorset
        self.colorset.reset()
        # init the rcParams
        if self.WHR > 1:
            plt.rcParams['figure.subplot.left']     = 0.1/self.WHR
            plt.rcParams['figure.subplot.right']    = 1 - 0.1/self.WHR
            plt.rcParams['figure.subplot.top']      = 0.9
            plt.rcParams['figure.subplot.bottom']   = 0.1
        else:
            plt.rcParams['figure.subplot.left']     = 0.1
            plt.rcParams['figure.subplot.right']    = 0.9
            plt.rcParams['figure.subplot.top']      = 1 - 0.1 * self.WHR
            plt.rcParams['figure.subplot.bottom']   = 0.1 * self.WHR
        # figure and axes
        self.figure:Figure = plt.figure(figsize=figsize)
        # if only 1 SumStats, or use circular layout
        if self.multitest == 1 or self.layout in ['cavern','firework'] or self.style in ['symmetric','overlap']:
            self.ax:list = [self.figure.subplots(1,1)]
        # if mulit SumStats
        else:
            if self.layout == 'vertical':
                self.ax:list        = self.figure.subplots(1, self.multitest)
                self.xzoom:float    = self.xzoom / self.multitest 
            else:
                self.ax:list        = self.figure.subplots(self.multitest, 1)
                self.yzoom:float    = self.yzoom / self.multitest
        if self.style == 'classic':
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
                if self.layout == 'horizontal':
                    self.plot_maht_horiz(axes=self.ax[idx],sumstats=plot_sumstats,quick=quick)
                elif self.layout == 'vertical':
                    self.plot_maht_vert(axes=self.ax[idx],sumstats=plot_sumstats,quick=quick)
        elif self.style == 'overlap':
            # axes to plot
            axex = self.ax[0]
            # get the max ymax
            ymax = -np.log10(self.threshold)
            for SS in self.sumstats:
                SS:SummaryStats
                # get ymax of all SS
                ymax_this = SS.data['log10_pvalue'].max() * 1.2
                if ymax_this > ymax: 
                    ymax = ymax_this
            # plot the axis
            # TODO: 兼容对称和叠加模式
            self.plot_classic_axis_horiz(axes=axex,sumstats=self.sumstats[0],ymax=ymax)
            for SS in self.sumstats:
                SS:SummaryStats
                # get ymax of all SS
                self.plot_classic_scattor_horiz(axes=axex,sumstats=SS,quick=1)
                if ymax_this > ymax: 
                    ymax = ymax_this
                
            # if self.
            # self.plot_maht_
        elif self.style == 'symmetric':
            pass
            # TODO 其他情况
        return self
    def plot_maht_vert(self):
        pass
    def plot_classic_axis_horiz(
            self,
            axes:Axes,
            sumstats:SummaryStats,
            ymax:float,
            ):
        # yunit
        yunit = ymax/150
        ymin = - 5 * yunit
        # ylim = ymax - ymin
        axes.set_ylim(ymin,ymax)
        # x axis:
        if sumstats.chrom > 0:
            x_min = sumstats.theta(...,sumstats.from_bp)
            x_max = sumstats.theta(...,sumstats.to_bp)
            x_sep = np.linspace(sumstats.from_bp,sumstats.to_bp,11)
            # -----
            axes.plot([x_min,x_max],[0,0],color='k',zorder=3)
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
                    verticalalignment   = 'top',
                    zorder = 3
                    )
            # CHROM
            axes.text(
                x   = np.pi,
                y   = -10*yunit / self.yzoom,
                s   = 'Chr%d(M)'%sumstats.chrom,
                fontsize = 14,
                horizontalalignment = 'center',
                verticalalignment   = 'center',
                zorder = 3
                )
        else:
            n_chrom:int = sumstats.data['chrom'].max()
            # simple line
            if self.style in ['overlap','symmetric']:
                axes.plot([x_min,x_max],[0,0],color='k',zorder=3)
                for chrom in range(1,n_chrom+1):
                    x_begin:float   = sumstats.theta(chrom, 0)
                    axes.plot([x_begin,x_begin],[2*yunit,-2*yunit],color='k')
            # color bar
            elif self.style == 'classic':
                self.colorset.reset()
                for chrom in range(1,n_chrom+1): #bar
                    x_begin:float   = sumstats.theta(chrom, 0)
                    x_end:float     = sumstats.theta(chrom, sumstats.chr_len[chrom])
                    x_width:float   = x_end - x_begin
                    x_loc:float     = (x_begin + x_end)/2
                    # color of this chromosome
                    color_this = self.colorset.next()
                    # chromosome range
                    axes.add_patch(
                        Rectangle(
                            xy      = (x_begin, -3 * yunit),
                            width   = x_width,
                            height  = 3 * yunit,
                            ec      = color_this[:3]/2,
                            color   = color_this,
                            zorder  = 3
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
                    # label of chrom on background
                    chrom_label = (lambda C:str(C) if C <= sumstats.heterosome['autosome'] else sumstats.heterosome[C])(chrom)
                    vert_horz = chrom<10 or chrom>sumstats.heterosome['autosome']
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
                        verticalalignment   = 'top',
                        zorder = 3
                        )
            # chromosome label
            for chrom in range(1,n_chrom+1):
                x_begin:float   = sumstats.theta(chrom, 0)
                x_end:float     = sumstats.theta(chrom, sumstats.chr_len[chrom])
                x_loc:float     = (x_begin + x_end)/2
                chrom_label = (lambda C:str(C) if C <= sumstats.heterosome['autosome'] else sumstats.heterosome[C])(chrom)
                axes.text(
                    x       = x_loc,
                    y       = - 4 * yunit,
                    s       = chrom_label,
                    fontsize= 14,
                    horizontalalignment = 'center',
                    verticalalignment   = 'top',
                    zorder = 3
                    )
            # CHROM
            axes.text(
                x   = np.pi,
                y   = -10*yunit/self.yzoom,
                s   = 'CHROM',
                fontsize = 14,
                horizontalalignment = 'center',
                verticalalignment   = 'center',
                zorder = 4
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
        axes.plot([-0.003*self.radian,-0.003*self.radian],[0,ymax],color='k',zorder = 3)
        # _|
        axes.plot([self.radian*1.003,self.radian*1.003],[0,ymax],color='k',zorder = 3)
        # _|
        axes.plot([self.radian*0.9985,self.radian*1.003],[0,0],color='k',zorder = 3)
        # --
        axes.plot([-0.003*self.radian,1.003*self.radian],[ymax,ymax],color='k',zorder = 3)
        axes.text(
            x   = (-0.022-0.01*np.ceil(np.log10(i))) * self.radian / self.xzoom,
            y   = ymax/2,
            s   = '-log10(P)',
            fontsize    = 14,
            rotation    = 90,
            horizontalalignment = 'center',
            verticalalignment   = 'center',
            zorder  = 4
            )
        axes.text(
            x   = self.radian/2,
            y   = ymax,
            s   = sumstats.name+': %d'%sumstats.data['size'].median(),
            fontsize    =16,
            horizontalalignment = 'center',
            verticalalignment   = 'bottom',
            zorder = 4
            )
        # significant threshold line
        axes.plot(
            np.linspace(0,self.radian,360),
            -np.log10(self.threshold)*np.ones(360),
            linestyle   = '-',
            c           = 'k',
            zorder      = 1,
            lw          = 0.5
            )
        axes.plot(
            np.linspace(0,self.radian,360),
            -np.log10(self.threshold/self.multitest)*np.ones(360),
            linestyle   = '--',
            c           = 'k',
            zorder      = 1,
            lw          = 0.5
            )
    def plot_classic_scattor_horiz(
            self,
            axes:Axes,
            sumstats:SummaryStats,
            quick:bool = True,
            color_mode:str = 'multi'
        ):
        not_significant = sumstats.significant(
            level       = 'negative',
            threshold   = self.threshold,
            mulitest    = self.multitest,
            window      = self.window
            )
        if quick:
            if len(not_significant) < 50000:
                not_significant = not_significant.copy()
            else:
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
        # get theta
        not_significant['theta']    = sumstats.theta(not_significant['chrom']   ,not_significant['pos']   )
        genome_significant['theta'] = sumstats.theta(genome_significant['chrom'],genome_significant['pos'])
        study_significant['theta']  = sumstats.theta(study_significant['chrom'] ,study_significant['pos'] )
        if sumstats.chrom > 0:
            axes.scatter(
                x       = not_significant['theta'],
                y       = not_significant['log10_pvalue'],
                s       = 20,
                alpha   = 1,
                color   = self.colorset.next(),
                zorder  = 2
                )
        else:
            n_chrom:int = sumstats.data['chrom'].max()
            self.colorset.reset()
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
                    color   = self.colorset.next(),
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
    def plot_maht_horiz(
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
        # theta
        sumstats.get_theta(chr_sep=self.chr_sep,radian=self.radian)
        # plot outframe
        self.plot_classic_axis_horiz(
            axes        = axes,
            sumstats    = sumstats,
            ymax        = ymax
        )
        # plot scatter
        if self.marker == 'scatter':
            self.plot_classic_scattor_horiz(
                axes        = axes,
                sumstats    = sumstats,
                quick       = quick
                )
        # gene text
        if self.showgene != 'False':
            self.mark_gene_horiz(
                axes        = axes,
                sumstats    = sumstats,
                ymax        = ymax
                )
    # mark up gene pos
    def mark_gene_horiz(self,axes:Axes,sumstats:SummaryStats,ymax:float):
        gene_to_show = sumstats.get_gene(
            level       = self.showgene,
            threshold   = self.threshold,
            window      = self.window
            ).sort_values('log10_pvalue').dropna(subset='gene')
        gene_to_show['theta'] = sumstats.theta(
            gene_to_show['chrom'],
            gene_to_show['pos']
        )
        gene_to_show['gene_y_loc'] = gene_to_show['log10_pvalue']
        i=0
        # the real yloc of gene text
        for key,values in gene_to_show.iterrows():
            text_yloc = gene_to_show['gene_y_loc'].values
            x_gene_protect_distance = self.radian * 0.04 * self.geneXdist / self.xzoom
            text_yloc[i] = gene_to_show[
                (gene_to_show['theta'] >= values['theta']  - x_gene_protect_distance)
                &(gene_to_show['theta'] <= values['theta'] + x_gene_protect_distance)
                ]['gene_y_loc'].max() + ymax * 0.025 / self.yzoom
        # Update the position of 'y' in the data table 
        # to ensure that the value of 'max' will change during the next loop, 
        # achieving the effect of shifting the label up.
            gene_to_show['gene_y_loc'] = text_yloc
            i += 1
            try:
                # if values.gene in self.:
                #     cr = 'k'
                # else:
                #     cr = 'r'
                axes.text(
                    x       = values['theta'],
                    y       = gene_to_show.loc[key,'gene_y_loc'],
                    s       = values.gene,
                    color   = 'r',
                    fontdict={
                        'fontstyle': 'italic',
                        'size'     : 13
                        },
                    horizontalalignment ='center',
                    verticalalignment   ='center'
                    )
            except:
                continue
    # get the rotation angle of chromosome text label
    def rotate_text(self,theta):
        if theta <= np.pi/self.radian :
            return theta
        else:
            return theta - np.pi/self.radian
    # save fig
    def save(self,path,**option):
        self.figure.savefig(fname=path,**option)