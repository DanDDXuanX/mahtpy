#!/usr/bin/env python3
# coding: utf-8
# lilinxuan@genomics.cn
# draw Manhattan plot py matplotlib

# using: python mahtplot_classic.py -i ifile.linear -o ofile.png
# suggestion: ifile can use the raw output of "plink --glm", 
# however using "grep 'ADD' ifile" to exclude the covariants rows can significantly reduce running time.

# v1.1 add axis name title ,debug cut10 function
# v1.2 add "TEST = ADD" 
# v1.3 allow argparse
# v1.4 gene name fontsize
# v1.5 scatter xdistance
# v1.6 gene?

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sys
import argparse
#import os
#import time
# 全局变量

cmap_choices=['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 
               'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 
               'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 
               'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 
               'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 
               'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 
               'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 
               'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 
               'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 
               'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 
               'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 
               'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 
               'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 
               'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 
               'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 
               'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 
               'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 
               'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 
               'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 
               'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 
               'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 
               'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 
               'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 
               'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 
               'viridis', 'viridis_r', 'winter', 'winter_r']

# 读取参数
parser = argparse.ArgumentParser(description="Draw Manhattan plot py matplotlib! v1.4 2021-12-06")
parser.add_argument("-i", "--input", help='The output file of plink --glm, .linear or .logistic',type=str)
parser.add_argument("-o", "--output", help='The output plot path, format is .png',type=str)
parser.add_argument("-t", "--threshold", help="Genome-wide significance threshold.", default=5e-8, type=np.float64)
#parser.add_argument("-c", "--cmap", help="The color schemes of mahtplot, these cmap is available: "+' '.join(cmap_choices), default='hsv', type=str)
parser.add_argument("-p", "--posinitial", help="The color pos of chr1 in plot.", default=0, type=np.int64) #v1.4
parser.add_argument("-d", "--distance", help="The distance between adjacent colors.", default=1, type=np.int64)
parser.add_argument("-n", "--numcolors", help="The number of colors in plot.", default=2, type=np.int64)
parser.add_argument("-s", "--showname", help="The name of traits on title.", default='', type=str)
# 多重比较的阈值
parser.add_argument("--multipletesting", help="Correct the threshold after multiple testing, by the times of GWAS tests.",default=1,type=int)
# col specify
parser.add_argument("--samplesize", help="Specify a columns of input file as samplesize of GWAS summary.",default='OBS_CT',type=str)
parser.add_argument("--chrom", help="Specify a columns of input file as chromosome index of GWAS summary.",default='#CHROM',type=str)
parser.add_argument("--position", help="Specify a columns of input file as base position of GWAS summary.",default='POS',type=str)
parser.add_argument("--pvalue", help="Specify a columns of input file as p-values of GWAS summary.",default='P',type=str)
#v1.5
parser.add_argument("-k", "--known", help="The exist known significant gene list of the trait", default='None', type=str)
parser.add_argument("--gap", help="The distance between chromosomes in plot, unit is bp.", default=20000000, type=np.int64)
parser.add_argument("--xProtDist", help="X-axis protection distance multiplier", default=1.0, type=np.float64)
parser.add_argument("--hidegene",help="Show annotate genes or not",action='store_true')

args = parser.parse_args()

#if args.cmap not in cmap_choices:
#    raise ValueError('Unrecognizable cmap name!'+args.cmap)

# col specify
SIZE = args.samplesize
PVAL = args.pvalue
CHR = args.chrom
POS = args.position

# 多重比较
mulitple = args.multipletesting

if type(args.input)!=str or type(args.output)!=str:
    raise ValueError('Miss ifile or ofile path!')

if (args.numcolors<0) or (args.distance<0) or (args.gap<0):
    raise ValueError('Invalid integer values!')

Dist = args.distance #临近配色的色彩距离，取1则是渐变颜色
Group = args.numcolors #绘制颜色的总数量
Init = args.posinitial #绘制的第一个颜色的位置
Chr_sep = args.gap #染色体之间的间隔距离
showname = args.showname # 展示在标题上的名字
#C_map = args.cmap #配色表
Threshold = args.threshold #显著性阈值
xProtDist = args.xProtDist #x轴保护距离乘数
hidegene = args.hidegene #是否显示基因标注

Dim = 2 * np.pi #在圆周上绘制的弧度#在classic里不建议修改
# 文件参数

ifile = args.input #输入的文件的路径
ofile = args.output #输入的图片的路径
known = args.known #已知基因名的路径

if known == 'None':
    mask_list = []
else:
    mask_list = open(known,'r').read().split('\n')#v1.5 新增一个masklist，用于给现存基因名标黑

if showname == '':
    showname = ofile.split('/')[-1].split('.')[0]

# 函数定义

# 处理X染色体数值位23
def X_23(x):
    if type(x) is int:
        return x
    else: # chr
        if x[0:3] == 'chr':
            x = x[3:]
        if x == "X":
            return 23
        else:
            return int(x) #1.5 返回值设置为int，混合的type导致（invalid ==）
X_23_uf=np.frompyfunc(X_23,1,1)

# 反转颜色
def reverse_c(x): #反转颜色
    global Dist,Group
    if type(x)!= int:
        try:
            x=int(x)
        except:
            raise TypeError('not supported datatype')
    return (Init+(x-1)%Group*Dist)%23
reverse_c_u=np.frompyfunc(reverse_c,1,1)

# 得到基因位置
def get_gene(chrom,pos): 
    if chrom == 23:
        chrom = 'chrX'
    else:
        chrom = 'chr'+str(chrom)
    Chr_b=gene_pos['#hg38.knownCanonical.chrom']==chrom
    Bg_b=gene_pos['hg38.knownCanonical.chromStart']<=pos
    Ed_b=gene_pos['hg38.knownCanonical.chromEnd']>=pos
    try:
        return gene_pos[Chr_b&Bg_b&Ed_b]['hg38.kgXref.geneSymbol'].values[0]
    except:
        return np.nan
get_gene_u=np.frompyfunc(get_gene,2,1)

# 达到180度之后翻转文本180度，更容易读,x_loc是在Dim中的比例
def rotate_text(x_loc): 
    if x_loc<= np.pi/Dim :
        return x_loc
    else:
        return x_loc - np.pi/Dim # np.pi/Dim 是图中展示的180度的位置

#载入指定色图、chr_len、gene_pos
#color_map=eval("plt.cm.{0}(np.linspace(0, 1, 23))".format(C_map))
data = np.array([(180/255, 216/255, 161/255,1,169/255, 160/255, 215/255,1)]).repeat(12,axis=0)
color_map = data.reshape(24,4)[0:23]
chr_len=pd.read_csv('./CHR.len',sep=' ',names=['chrom','len'],index_col='chrom')['len']
gene_pos=pd.read_csv('./KnownCanonicalGene.Drop.hg38.txt',sep='\t')

#载入输入文件
isSevere_logist=pd.read_csv(ifile,sep='\t',dtype={CHR:str})#v1.5 指定CHROM以字符串形式读取
sample_size=isSevere_logist[SIZE].median()#v1.4

#处理pvalue=0 v1.5
if isSevere_logist[PVAL].min() == 0:
    isSevere_logist[PVAL] = isSevere_logist[PVAL].replace({0:1e-300})

isSevere_logist[CHR]=X_23_uf(isSevere_logist[CHR])#处理掉X染色体
isSevere_logist = isSevere_logist[[CHR,POS,PVAL]].copy() #v1.2 添加了只筛选到TEST = ADD

chr_len_pre=chr_len[:-2].cumsum()+np.linspace(1,23,23,dtype=int)*Chr_sep
chr_len_pre.index=range(1,24)
chr_len_pre[0]=0
chr_len_total=chr_len_pre[23]#环上的染色体总长

#color_space=reverse_c_u(isSevere_logist['#CHROM']).astype(np.int64) #v1.5
color_index=reverse_c_u(np.linspace(1,23,23))

isSevere_logist['global_pos']=np.frompyfunc((lambda chrom,pos:pos+chr_len_pre[chrom-1]),2,1)(isSevere_logist[CHR].astype(np.int64),isSevere_logist[POS])
isSevere_logist['theta']=isSevere_logist['global_pos']/chr_len_total*Dim #环上的角度RAD
isSevere_logist['p_log_10']=-np.log10(isSevere_logist[PVAL])

#isSevere_logist['p_log_10']=isSevere_logist['p_log_10'].map(max_log10p_11) # cut到 max11

isSevere_logist['gene']=get_gene_u(isSevere_logist[isSevere_logist[PVAL]<=Threshold][CHR],isSevere_logist[isSevere_logist[PVAL]<=Threshold][POS])
not_significant = isSevere_logist[isSevere_logist[PVAL]>(Threshold)]
genome_significant=isSevere_logist[(isSevere_logist[PVAL]<=Threshold) & (isSevere_logist[PVAL]>(Threshold/8))]
study_significant=isSevere_logist[isSevere_logist[PVAL]<=(Threshold/8)]
significant = pd.concat([genome_significant,study_significant])
top_snps=pd.merge(significant.groupby('gene',as_index=False).max()[['gene','p_log_10']],significant,left_on=['gene','p_log_10'],right_on=['gene','p_log_10'],how='left').drop_duplicates(['gene','p_log_10'],keep='first').sort_values('p_log_10')
top_snps['text_y_loc']=top_snps['p_log_10']

# 绘图
# 参数
plt.rcParams['figure.subplot.left']=0.05
plt.rcParams['figure.subplot.right']=0.95
plt.rcParams['figure.subplot.top']=0.9
plt.rcParams['figure.subplot.bottom']=0.1

fig=plt.figure(figsize=(16,8))
ax=plt.subplot(111)

ax.axis('off')#隐藏坐标轴

# 设置y轴范围 #v1.4 上限+0.1x
ax.set_ylim(-isSevere_logist['p_log_10'].max()*0.03,max(-np.log10(Threshold),isSevere_logist['p_log_10'].max())*1.2)

#画曼哈顿图 #z_order 设置图层顺序
#maht=ax.scatter(x=isSevere_logist['theta'],y=isSevere_logist['p_log_10'],
#               s=10,alpha=0.5,
#               #c=color_space,cmap=C_map, #v1.4debug 
#               color=color_map[color_space], #这里有bug
#               zorder=2)
#v1.5 尝试分染色体逐个做散点
for i in range(0,23):
    point_this = not_significant[not_significant[CHR]==(i+1)]
    ax.scatter(x=point_this['theta'],y=point_this['p_log_10'],
               s=20,alpha=0.5,
               #c=color_space,cmap=C_map, #v1.4debug 
               color=color_map[color_index[i]], #这里有bug
               zorder=2)

# 显著的 标为红色
ax.scatter(x=genome_significant['theta'],y=genome_significant['p_log_10'],
           s=20,alpha=0.5,color=[0.2,1,0.2,1],zorder=2)
ax.scatter(x=study_significant['theta'],y=study_significant['p_log_10'],
           s=40,alpha=0.5,color=[1,0.2,0.2,1],zorder=2)


#显著相关阈值线
dash_line=ax.plot(np.linspace(0,Dim,360),-np.log10(Threshold)*np.ones(360),
                  linestyle='-',c='k',zorder=1,lw=0.5)
dash_line=ax.plot(np.linspace(0,Dim,360),-np.log10(Threshold/mulitple)*np.ones(360),
                  linestyle='--',c='k',zorder=1,lw=0.5)

ylim = ax.get_ylim()[1]-ax.get_ylim()[0]
ymax = ax.get_ylim()[1]

# x轴
for i in range(0,23): #bar
    x_loc=(chr_len_pre[i]+chr_len[i]/2)/chr_len_total
    chrom_label_rotate=rotate_text(x_loc) 
    ax.bar(
        x=x_loc*Dim,
        height=ylim*3/150,#色条高度 #v1.4 色条的宽度削减1/3
        width=chr_len[i]/chr_len_total*Dim,
        bottom = -ylim/50,ec=color_map[color_index[i]][:3]/2,
        color=color_map[color_index[i]],
        zorder=3
           )
    ax.bar(
        x=x_loc*Dim,
        height=ymax,
        width=chr_len[i]/chr_len_total*Dim,
        bottom = 0,
        color=color_map[color_index[i]],
        zorder=3,
        alpha=0.1
        )
    ax.text(
        x=(chr_len_pre[i]+chr_len[i]/2)/chr_len_total*Dim,
        y=-ylim*4/150,
        s=chr_len.index[i][3:],
        fontsize=14,
        horizontalalignment='center',verticalalignment='top'
        )#染色体名称
    if i < 9 or i == 22:
        ax.text(
            x=(chr_len_pre[i]+chr_len[i]/2)/chr_len_total*Dim,
            y=ymax-ylim*4/150,
            s=chr_len.index[i][3:],
            fontdict={'size':24,'weight':'bold','color':color_map[color_index[i]]},
            color=color_map[color_index[i]],
            alpha=0.2,
            horizontalalignment='center',verticalalignment='top'
            )#染色体名称
    else:
        ax.text(
            x=(chr_len_pre[i]+chr_len[i]/2)/chr_len_total*Dim,
            y=ymax-ylim*4/150,
            s=chr_len.index[i][3:],
            fontdict={'size':16,'weight':'bold','color':color_map[color_index[i]]},
            color=color_map[color_index[i]],
            alpha=0.2,
            rotation=90,
            horizontalalignment='center',verticalalignment='top'
            )#染色体名称
ax.text(x=np.pi,y=-ylim/15,s='CHROM',fontsize=14,horizontalalignment='center',verticalalignment='center')

#gene text fontsize：8-x保护范围0.2-y保护范围0.02x; v1.4:fz12-x保护范围0.3-y保护范围0.027x(贪婪至13,x贪婪至0.2)
if hidegene == False:
    i=0
    for key,values in top_snps.iterrows():
        tyl=top_snps.text_y_loc.values
        tyl[i]=top_snps[(top_snps.theta>=values.theta-0.2*xProtDist)&(top_snps.theta<=values.theta+0.2*xProtDist)].text_y_loc.max()+ylim*0.025
        top_snps.text_y_loc=tyl #更新数据表中的y的位置，保证下一次循环的时候的max的值会被改变，从而达到标签上移的效果
        i+=1
        try:
            if values.gene in mask_list:
                cr = 'k'
            else:
                cr = 'r'#v1.5 在masklist里面的基因名标记为黑色，反之为红色
            ax.text(x=values.theta,y=top_snps.loc[key].text_y_loc,
                    s=values.gene,color=cr,
                    fontsize=13,
                    horizontalalignment='center',verticalalignment='center')
        except:
            continue
else:
    pass

# y轴
for g in range(1,100):
    if int(ymax/g)<=10:
        break
    else:
        pass

# =|_
for i in range(0,int(ymax+1),g):
    ax.text(x=-0.06,y=i,s=i,fontsize=14,
            horizontalalignment='right',verticalalignment='center')#y轴坐标轴值
    ax.plot([-0.05,-0.02],[i,i],color='k')
# |_
ax.plot([-0.02,-0.02],[0,ymax],color='k')
# _|
ax.plot([Dim+0.02,Dim+0.02],[0,ymax],color='k')
# _|
ax.plot([Dim-0.01,Dim+0.02],[0,0],color='k')
# --
ax.plot([-0.02,Dim+0.02],[ymax,ymax],color='k')

ax.text(x=-0.14-0.06*np.ceil(np.log10(i)),y=ymax/2,s='-log10(P)',
        fontsize=14,horizontalalignment='center',verticalalignment='center',rotation=90)
ax.text(x=np.pi,y=ymax,s=showname+': '+('%d'%sample_size),
        fontsize=16,horizontalalignment='center',verticalalignment='bottom')#v1.4 在title上标注samplesize，fz+2

fig.savefig(ofile,format='png',dpi=150) #保存
