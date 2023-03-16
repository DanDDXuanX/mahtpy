#!/usr/bin/env python3
# coding: utf-8
# lilinxuan@genomics.cn
# summary statistics class definition

import pandas as pd
import numpy as np
from pathlib import Path

class SumstatsError(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class SummaryStats:
    # class variables
    chr_len_dict:dict = {
        'hg38' : Path('./assets/CHR.len')
    }
    gene_pos_dict:dict = {
        'hg38' : Path('./assets/KnownCanonicalGene.Drop.hg38.txt')
    }
    # init with a file, pd.DF or a np.2darray
    def __init__(
            self,
            file_input  : str           = None,
            frame_input : pd.DataFrame  = None,
            array_input : np.ndarray    = None,
            seperator   : str           = '\t',
            col_specify : dict          = {},
            ref_genome  : str           = 'hg38'
            ) -> None:
        # load col specify
        col_name = {
                'chrom'     : '#CHROM',
                'pos'     : 'POS',
                'pvalue'  : 'P',
                'size'    : 'OBS_CT' 
                }
        for key in col_specify.keys():
            if key not in col_name.keys():
                continue
            else:
                col_name[key] = col_specify[key]
        # load summary from file
        if file_input is not None:
            filepath = Path(file_input)
            try:
                # an excel table file is specified
                if filepath.suffix in ['xlsx','xls']:
                    self.data = pd.read_excel(filepath)
                # load as csv
                else:
                    self.data = pd.read_csv(filepath,sep=seperator)
            except Exception as E:
                raise SumstatsError(
                    "Unable to load file '{}' as a DataFrame, due to: {}"
                    .format(filepath,E)
                )
            # n col is 1, bad seperator.
            if self.data.shape[1] == 1:
                raise SumstatsError(
                    "Error occurred while seperating columns with seperator '{}'"
                    .format(seperator)
                )
        # load summary from DF
        elif frame_input is not None:
            self.data = frame_input
        # load summary from array
        elif array_input is not None:
            # check if array is 2 dim
            if array_input.ndim != 2:
                raise SumstatsError(
                        "Dimension of array like input is not 2."
                    )
            else:
                self.data = pd.DataFrame(array_input)
        # No input Error
        else:
            raise SumstatsError(
                "Need an input to initialize SummaryStats."
            )
        # keep the specified cols
        try:
            self.data = self.data[list(col_name.values())].copy()
        except KeyError as E:
            raise SumstatsError(
                "Specified columns: {}".format(E)
            )
        # change the col name of table
        self.data.rename(
            columns = pd.Series(data=col_name.keys(),index=col_name.values()),
            inplace = True
        )
        # convert the dtype and format of 'chrom' col
        self.data['chrom'] = self.reformat_chrom(self.data['chrom'])
        # remove the invalid value
        try:
            self.data = self.data[self.data['chrom'] != -1].copy()
            self.data['pos'] = self.data['pos'].astype(int)
            self.data['size'] = self.data['size'].astype(int)
            self.data['pvalue'] = self.data['pvalue'].astype(float).replace({0:1e-300})
        except ValueError as E:
            SumstatsError(
                "Summary Table contain invalid value, due to: {}"
                .format(E)
            )
        # calculate the log 10 pvalue
        self.data['log10_pvalue'] = -np.log10(self.data['pvalue'])
        # load chrlen and gene pos of specified ref genome
        if ref_genome not in self.chr_len_dict.keys():
            raise SumstatsError(
                "The specified reference genome {} is not available."
                .format(ref_genome)
            )
        else:
            self.ref_genome = ref_genome
            self.load_info()
    # load chrlen and gene pos of specified ref genome
    def load_info(self)->None:
        # length of chrom
        if type(self.chr_len_dict[self.ref_genome]) is not pd.Series:
            readfile = pd.read_csv(
                self.chr_len_dict[self.ref_genome],
                sep   = ' ',
                names = ['chrom','len'],
            )
            readfile['chrom'] = self.reformat_chrom(readfile['chrom'])
            self.chr_len_dict[self.ref_genome] = readfile.set_index('chrom')['len']
        else:
            pass
        self.chr_len:pd.Series = self.chr_len_dict[self.ref_genome]
        # position of known gene on ref genome
        if type(self.gene_pos_dict[self.ref_genome]) is not pd.DataFrame:
            readfile = pd.read_csv(
                self.gene_pos_dict[self.ref_genome],
                sep = '\t'
            )
            readfile['#hg38.knownCanonical.chrom'] = self.reformat_chrom(
                readfile['#hg38.knownCanonical.chrom']
            )
            self.gene_pos_dict[self.ref_genome] = readfile
        else:
            pass
        self.known_gene:pd.DataFrame = self.gene_pos_dict[self.ref_genome]
    # calculate the display position of variants
    def get_theta(
            self,
            chr_sep:int     = 20000000, 
            radian:float    = 2*np.pi, 
            chrom:int       = -1,
            from_bp:int     = -1,
            to_bp:int       = -1,
            locus:int       = -1,
            threshold       = 5e-8,
            window          = 500000
        )->pd.Series:
        """
        calculate the display position in mathplot of each variants.\n

        Parameters:
        ----------
            chr_sep : int
                distance between chromosomes, unit is bp, default value is 20M.
            radian  : float
                total radian, or range in x axis to plot mahtplot, default value is 2*pi
            chrom   : int
                only display the specified chromosome, default is -1, means show all chromosomes.
            from_bp : int
            to_bp   : int
                only effective with chrom specified. the start and end pos to display.
            locus   : int
                only display the specified locusID, default is -1, means inoperative.
                Note: specified locus has higher priority than specified chromosome
            threshold : float
                genome-wide significant threshold, default is 5e-8.
            window      : int
                significant SNPs which distance less than window are merged into a loci, window default is 500kb
        Returns:
        ----------
            DataFrame
                SumStats with theta (postion to plot) value of variants in specified range.
        """
        # specify a locus to display
        if locus != -1:
            locus_top   = self.significant(threshold=threshold)
            locus_top   = locus_top.loc[locus_top.query("locusID==@locus")['log10_pvalue'].idxmax()]
            # return values
            chrom   = locus_top['chrom']
            from_bp = locus_top['pos'] - window
            to_bp   = locus_top['pos'] + window
            # sumstats range to display
            sumstats_to_display = (
                self.data[self.data['chrom']==chrom]
                .query("pos >= @from_bp and pos < @to_bp")
            )
            chr_len_total = 2 * window
            theta:pd.Series = (
                (sumstats_to_display['pos'] - from_bp)
                / chr_len_total * radian
            )
        # specify a chrom range to display
        elif chrom != -1:
            if from_bp == -1:
                from_bp = 0
            if to_bp == -1:
                to_bp = self.chr_len[chrom]
            sumstats_to_display = (
                self.data[self.data['chrom']==chrom]
                .query("pos >= @from_bp and pos < @to_bp")
            )
            chr_len_total = to_bp - from_bp
            theta:pd.Series = (
                (sumstats_to_display['pos'] - from_bp)
                / chr_len_total * radian
            )
        # display whole genome
        else:
            sumstats_to_display = self.data.copy()
            n_chrom:int = sumstats_to_display.max()
            # number of bp before each chrom
            chr_len_pre:pd.Series = (
                self.chr_len[:n_chrom].cumsum() + 
                np.arange(1,n_chrom+1,1,dtype=int) * chr_sep
                )
            chr_len_pre.index = np.arange(1,n_chrom+1,1,dtype=int)
            chr_len_pre[0]    = 0
            chr_len_total = chr_len_pre[n_chrom]
            # get the global pos (theta) of each variants
            theta:pd.Series = (
                np.frompyfunc((lambda chrom,pos : pos + chr_len_pre[chrom-1]), 2, 1)
                (
                    sumstats_to_display['chrom'],
                    sumstats_to_display['pos']
                )
                / chr_len_total * radian
            )
        # return
        sumstats_to_display['theta'] = theta
        return sumstats_to_display, [chrom,from_bp,to_bp]
    # convert chrom col as int type
    def reformat_chrom(self,chr_col:pd.Series)->pd.Series:
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
    # get the genome wide significant, study wide significant, and not significant variants.
    def significant(
            self,
            level:str       = 'study',
            threshold:float = 5e-8,
            mulitest:int    = 1,
            window:int      = 500000         
            ) -> pd.DataFrame:
        """
        get a subset of summary variants,
        at threshold of genome wide significant, 
        study wide significant, or not significant.

        Parameters:
        ----------
            level : str
                level of significance, possible value in ['negative','genome','study'], default is 'study'
                Note:   level 'genome' will return these variants pass genome-wide threshold, 
                        but not pass study-wide threshold!
            threshold : float
                genome-wide significant threshold, default is 5e-8.
            mulitest : int
                mulitple test time, to adjust study-wide threshold, default is 1.
            window      : int
                significant SNPs which distance less than window are merged into a loci, window default is 500kb
        Returns:
        ----------
            DataFrame
                subset of variants , at specified significant level.
        """
        # if get negative
        if level == 'negative':
            return self.data.query("pvalue > @threshold")
        # multiple testing threshold
        mulitest_threshold = threshold / mulitest
        # get locusIDS
        significants:pd.DataFrame = (
            self.data
            .query( "pvalue <= @threshold" )
            .sort_values(by=['chrom','pos'])
            .copy()
        )
        locusID = 0
        last = None
        for key,this in significants.iterrows():
            if last is None:
                pass
            else:
                if last['chrom']!=this['chrom']:
                    locusID += 1
                elif this['pos'] - last['pos'] > window:
                    locusID += 1
                else:
                    pass
            significants.loc[key,'locusID'] = locusID
            last = this
        self.window = window
        significants['locusID'] = significants['locusID'].astype(int)
        # return significant variants
        if level == 'genome':
            return significants.query("pvalue > @mulitest_threshold")
        elif level == 'study':
            return significants.query("pvalue <= @mulitest_threshold")
        else:
            raise SumstatsError("Invalid significant level: '{}'".format(level))
    # get mapped gene of variants.
    def get_gene(
            self,
            threshold:float = 5e-8,
            level:str = 'snp',
            window:int = 500000,
            ) -> pd.Series:
        """
        get the mapped gene names of each significant variant,

        Parameters:
        ----------
            threshold   : float
                genome-wide significant threshold, default is 5e-8.
            level       : str
                possible value in ['snp','gene','loci'],\n
                if level is 'snp', return all significant SNPs with mapped genes annotated,\n
                if level is 'gene', only return the top SNPs of each mapped gene,\n
                if level is 'loci', adjacent significant SNPs are merged into a loci, only return the top snps of each loci.\n
                default is 'snp'.
            window      : int
                significant SNPs which distance less than window are merged into a loci, window default is 500kb
        
        Returns:
        ----------
            Series
                mapped gene name list of variants which reached the threshold.
        """
        # get mapped gene from 
        def mapped(chrom,pos): 
            # if chrom == 23:
            #     chrom = 'chrX'
            # else:
            #     chrom = 'chr'+str(chrom)
            Chr_b=self.known_gene['#hg38.knownCanonical.chrom']==chrom
            this_chrom = self.known_gene[Chr_b]
            Bg_b = this_chrom['hg38.knownCanonical.chromStart']<=pos
            Ed_b = this_chrom['hg38.knownCanonical.chromEnd']>=pos
            try:
                return this_chrom[Bg_b&Ed_b]['hg38.kgXref.geneSymbol'].values[0]
            except:
                return np.nan
        # get all significant snps
        sig = self.significant(threshold=threshold,window=window).copy()
        if level == 'loci':
            sig = sig.loc[sig.groupby('locusID')['log10_pvalue'].idxmax()].copy()
        # get all gene annotats
        sig['gene'] = (
            np.frompyfunc(mapped, 2, 1)(
                sig['chrom'], sig['pos']
            )
        )
        # case : level
        if level == 'snp':
            return sig
        elif level == 'gene':
            return (
                sig
                .dropna(subset='gene')
                .loc[
                    sig
                    .groupby('gene')['log10_pvalue']
                    .idxmax()
                ] 
            )
        elif level == 'loci':
            return sig
        else:
            raise SumstatsError("Invalid gene annotate level: '{}'".format(level))
