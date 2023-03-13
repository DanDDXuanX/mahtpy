#!/usr/bin/env python3
# coding: utf-8
# lilinxuan@genomics.cn
# summary statistics class definition

import pandas as pd
import numpy as np
from pathlib import Path

class MahtplotError(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class SummaryStats:
    # class variables
    chr_len_dict:dict = {
        'hg38' : Path('./CHR.len')
    }
    gene_pos_dict:dict = {
        'hg38' : Path('./KnownCanonicalGene.Drop.hg38.txt')
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
                raise MahtplotError(
                    "Unable to load file '{}' as a DataFrame, due to: {}"
                    .format(filepath,E)
                )
            # n col is 1, bad seperator.
            if self.data.shape[1] == 1:
                raise MahtplotError(
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
                raise MahtplotError(
                        "Dimension of array like input is not 2."
                    )
            else:
                self.data = pd.DataFrame(array_input)
        # No input Error
        else:
            raise MahtplotError(
                "Need an input to initialize SummaryStats."
            )
        # keep the specified cols
        try:
            self.data = self.data[list(col_name.values())].copy()
        except KeyError as E:
            raise MahtplotError(
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
            MahtplotError(
                "Summary Table contain invalid value, due to: {}"
                .format(E)
            )
        # calculate the log 10 pvalue
        self.data['log10_pvalue'] = -np.log10(self.data['pvalue'])
        # load chrlen and gene pos of specified ref genome
        if ref_genome not in self.chr_len_dict.keys():
            raise MahtplotError(
                "The specified reference genome {} is not available."
                .format(ref_genome)
            )
        else:
            self.ref_genome = ref_genome
            self.load_info()
    # load chrlen and gene pos of specified ref genome
    def load_info(self)->None:
        # length of chrom
        if self.chr_len_dict[self.ref_genome] is not pd.Series:
            self.chr_len_dict[self.ref_genome] = pd.read_csv(
            self.chr_len_dict[self.ref_genome],
            sep       = ' ',
            names     = ['chrom','len'],
            index_col = 'chrom'
            )['len']
        else:
            pass
        self.chr_len:pd.Series = self.chr_len_dict[self.ref_genome]
        # position of known gene on ref genome
        if self.gene_pos_dict[self.ref_genome] is not pd.DataFrame:
            self.gene_pos_dict[self.ref_genome] = pd.read_csv(
            self.gene_pos_dict[self.ref_genome],
            sep = '\t'
            )
        else:
            pass
        self.known_gene:pd.DataFrame = self.gene_pos_dict[self.ref_genome]
    # calculate the display position of variants
    def get_theta(self, chr_sep:int=20000000, radian:float=2*np.pi)->pd.Series:
        """
        calculate the display position in mathplot of each variants.\n

        Parameters:
        ----------
            chr_sep : int
                distance between chromosomes, unit is bp, default value is 20M.
            radian  : float
                total radian, or range in x axis to plot mahtplot, default value is 2*pi
        
        Returns:
        ----------
            Series
                theta (postion to plot) value of each variants.
        """
        n_chrom:int = self.data['chrom'].max()
        # number of bp before each chrom
        chr_len_pre:pd.Series = (
            self.chr_len[:n_chrom].cumsum() + 
            np.arange(1,n_chrom+1,1,dtype=int) * chr_sep
            )
        chr_len_pre.index = np.arange(1,n_chrom+1,1,dtype=int)
        chr_len_pre[0]    = 0
        chr_len_total = chr_len_pre[n_chrom]
        # get the global pos (theta) of each variants
        self.data['theta'] = (
            np.frompyfunc((lambda chrom,pos : pos + chr_len_pre[chrom-1]), 2, 1)
            (
                self.data['chrom'],
                self.data['pos']
            )
            / chr_len_total * radian
        )
        # return
        return self.data['theta']
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
                    else:
                        return int(x)
                # not supported type
                else:
                    return -1
            except ValueError:
                return -1
        return np.frompyfunc(reformat,1,1)(chr_col)
    # get the genome wide significant, study wide significant, and not significant variants.
    def significant(
            self,
            level:str       = 'study',
            threshold:float = 5e-8,
            mulitest:int    = 1
            ) -> pd.DataFrame:
        """
        get a subset of summary variants,
        at threshold of genome wide significant, 
        study wide significant, and not significant.

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
        
        Returns:
        ----------
            DataFrame
                subset of variants, at specified significant level.
        """
        mulitest_threshold = threshold / mulitest
        if level == 'negative':
            return self.data.query("pvalue > @threshold")
        elif level == 'genome':
            return self.data.query("pvalue <= @threshold and pvalue > @mulitest_threshold")
        elif level == 'study':
            return self.data.query("pvalue <= @mulitest_threshold")
        else:
            raise MahtplotError("Invalid significant level: '{}'".format(level))