#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 18:04:47 2018

@author: valeria
"""

import numpy as np
import math
import os
import gdal 
from gdalconst import * 
from osgeo import osr
import matplotlib
matplotlib.use('qt5agg')
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import datetime
import sys
sys.path.append('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')

#PLOTTING
from matplotlib import colors
class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)
    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

import Read_Hi_data


def find_absdif(Phi,Chi):
    param = Phi-Chi
    return param

def plot_map_Greenland(param, clon, clat, outfname, title):
    m = Basemap(resolution="i",
                projection='laea', lat_ts=90, lat_0=90., lon_0=0.,
                llcrnrlon= clon[-1,0], llcrnrlat= clat[-1,0] ,
                urcrnrlon= clon[0,-1], urcrnrlat= clat[0,-1])

    param_nn = param.copy()
    param_nn[np.isnan(param)]=0

    cmap = 'RdBu_r'
    # m.imshow(param_nn, origin = 'upper', norm=MidpointNormalize(midpoint=0, vmax=2, vmin =-7), cmap = cmap)
    m.imshow(param_nn, origin='upper', vmin= -5, vmax = 5, cmap= discrete_cmap(10, base_cmap=cmap))
    m.drawcoastlines()
    plt.colorbar()
    m.drawmeridians(np.arange(-180,180,10))
    m.drawparallels(np.arange(55,80,10))
    plt.title(title)
    plt.savefig(outfname)
    plt.close()
    return

INDIR_Cr = '/home/lera/NIERSC/DATA/Cryosat'
INDIR_PIOMAS = '/home/lera/NIERSC/DATA/PIOMAS/EASE_grid'
flist_Cr =  [os.path.join(INDIR_Cr, f) for f in os.listdir(INDIR_Cr)]
flist_Cr.sort()
flist_PIOMAS = [os.path.join(INDIR_PIOMAS, f) for f in os.listdir(INDIR_PIOMAS)]
flist_PIOMAS.sort()
OUTDIR = '/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/figs/Greenland/'
abs_dif = []


for i in range(len(flist_PIOMAS)):
    print i
    date_PIOMAS = Read_Hi_data.extract_date_PIOMAS(flist_PIOMAS[i])
    for j in range(len(flist_Cr)):
        date_Cr = Read_Hi_data.extract_date_Cr(flist_Cr[j])
        if date_Cr==date_PIOMAS:
            clat, clon, chi, chi_u = Read_Hi_data.read_Cryosat(flist_Cr[j])
            plat, plon, phi = Read_Hi_data.read_PIOMAS(flist_PIOMAS[i])
            test = np.nanmin(phi)
            # if test<0:
            #     print date_PIOMAS
            data = phi-chi
            data = Read_Hi_data.extract_Greenland(data)
            data = Read_Hi_data.cut_domain(data)
            clat = Read_Hi_data.cut_domain(clat)
            clon = Read_Hi_data.cut_domain(clon)
            abs_dif.append(data)
            outfname = OUTDIR+str(date_PIOMAS.date())+'_abs_dif.png'
            title = str(date_PIOMAS.date())+'_abs_dif'
            plot_map_Greenland(data,clon,clat,outfname,title)
            

# v_min = []
# v_max = []
# std = []
# mean = []
#
# for i in range(len(abs_dif)):
#     v_min.append(np.nanmin(abs_dif[i]))
#     v_max.append(np.nanmax(abs_dif[i]))
#     std.append(np.nanstd(abs_dif[i]))
#     mean.append(np.nanmean(abs(abs_dif[i])))