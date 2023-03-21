#!/usr/bin/env python3
# coding: utf-8
# lilinxuan@genomics.cn

# available color map to plot.

from matplotlib.cm import _colormaps as mpl_cm
from matplotlib.colors import ColorConverter
import numpy as np

class ColorSetError(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class ColorSet:
    def __init__(
            self,
            colormap:str        =   None,
            colorset:list       =   [(0.706, 0.847, 0.631, 1.0),(0.663, 0.627, 0.843, 1.0)],
            **cmap_options:dict
            ) -> None:
        """
        Class ColorSet
        ----------
            a iterable object of color set.
            initialize by given a list of color, or name of matplotlib colormap.

        Parameters:
        ----------
            colorset : list
                list or array like, a list of color. default is [light-green, light-purple]\n
                each color element should be a valid matplotlib color value, for example:\n
                1. RGB(A) tuple : (1, 0.5, 0.5, 0.5)\n
                2. hex RGB string : "#0000ff"\n
                3. color name string : "red"\n
                more detail : https://matplotlib.org/stable/tutorials/colors/colors.html
            
            colormap : str
                a name of matplotlib colormap.\n
                more detail : https://matplotlib.org/stable/tutorials/colors/colormaps.html

            **cmap_option:
                keyword arguments to specify how to use the matplotlib cmap.

                ncolor : int
                    select number of color evenly from colormap.
                init : int
                    the first color index in selected colorset.
                dist : int
                    the index spacing from last output color.
            
        method:
        ----------
            next
                output a color from colorset.
        """
        # cmap options
        self.cmap_option = {
                'init'      :0,
                'dist'      :1,
                'ncolor'    :2, 
            }
        # colorset and colormap
        if colormap in mpl_cm.keys():
            for key in cmap_options:
                if key in self.cmap_option.keys():
                    self.cmap_option[key] = cmap_options[key]
                else:
                    pass
            cm = mpl_cm[colormap]
            self.colorset:np.ndarray = cm(np.linspace(0,1,self.cmap_option['ncolor']))
            self.this = self.cmap_option['init']
        elif len(colorset) != 0:
            color_converter = ColorConverter()
            self.colorset = []
            for C in colorset:
                try:
                    self.colorset.append(color_converter.to_rgba(c=C,alpha=1))
                except:
                    ColorSetError("Invalid color key '{}'.".format(C))
            self.colorset:np.ndarray = np.asarray(self.colorset)
            self.this = 0
            self.cmap_option = {
                'init'      :0,
                'dist'      :1,
                'ncolor'    :len(colorset), 
            }
        else:
            raise ColorSetError("Invalid arguments for colorset.")
    # get next color
    def next(self):
        # return in this iter
        to_return = self.colorset[self.this]
        # update this
        self.this = self.this + self.cmap_option['dist']
        # move this back in ncolor, to loop
        if self.this >= self.cmap_option['ncolor']:
            self.this = self.this % self.cmap_option['ncolor']
        return to_return
    # iterable
    def __iter__(self):
        return self
    def __next__(self):
        return self.next()
    # subscriptable
    def __getitem__(self,key:int):
        self.this = key
        return self.next()
