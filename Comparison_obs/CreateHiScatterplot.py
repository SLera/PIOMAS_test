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
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime
import sys

sys.path.append('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')

import Read_Hi_data


def extract_date_PIOMAS(path):
    return datetime.datetime.strptime(path[-9:-3], '%Y%m')


def extract_date_Cr(path):
    return datetime.datetime.strptime(path[-15:-9], '%Y%m')


def find_absdif(Phi, Chi):
    param = Phi - Chi
    return param

# def plot_map_Greenland(param, clon, clat, outfname, title):
#     m = Basemap(resolution="i",
#                 projection='laea', lat_ts=90, lat_0=90., lon_0=0.,
#                 llcrnrlon=clon[-1, 0], llcrnrlat=clat[-1, 0],
#                 urcrnrlon=clon[0, -1], urcrnrlat=clat[0, -1])
#
#     m.imshow(param, vmin=-10, vmax=2, origin='upper')
#     m.drawcoastlines()
#     plt.colorbar()
#     m.drawmeridians(np.arange(-180, 180, 10))
#     m.drawparallels(np.arange(55, 80, 10))
#     plt.title('PIOMAS-Cryosat' + title)
#     plt.savefig(outfname)
#     plt.close()
#     return


INDIR_Cr = '/home/valeria/DATA/Cryosat'
INDIR_PIOMAS = '/home/valeria/DATA/PIOMAS/v2.1/EASE_grid/rricker'
flist_Cr = [os.path.join(INDIR_Cr, f) for f in os.listdir(INDIR_Cr)]
flist_Cr.sort()
flist_PIOMAS = [os.path.join(INDIR_PIOMAS, f) for f in os.listdir(INDIR_PIOMAS)]
flist_PIOMAS.sort()
OUTDIR = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/figs/Greenland/'

months = np.array([10,11,12,1,2,3,4])

chi_months = []
phi_month = []

for i in range(len(flist_PIOMAS)):
    print i
    date_PIOMAS = extract_date_PIOMAS(flist_PIOMAS[i])
    if date_PIOMAS.month in months:
        for j in range(len(flist_Cr)):
            date_Cr = extract_date_Cr(flist_Cr[j])
            if date_Cr == date_PIOMAS:
                clat, clon, chi, chi_u = Read_Hi_data.read_Cryosat(flist_Cr[j])
                plat, plon, phi = Read_Hi_data.read_PIOMAS(flist_PIOMAS[i])
                ind_filter_P = np.isnan(chi)
                phi[ind_filter_P] = np.nan

                chi = Read_Hi_data.extract_Greenland(chi)
                chi_u = Read_Hi_data.extract_Greenland(chi_u)
                phi = Read_Hi_data.extract_Greenland(phi)

                phi = Read_Hi_data.cut_domain(phi)
                clat = Read_Hi_data.cut_domain(clat)
                clon = Read_Hi_data.cut_domain(clon)
                chi_u = Read_Hi_data.cut_domain(chi_u)
                chi = Read_Hi_data.cut_domain(chi)
                abs_dif.append(data)


v_min = []
v_max = []
std = []
mean = []

for i in range(len(abs_dif)):
    v_min.append(np.nanmin(abs_dif[i]))
    v_max.append(np.nanmax(abs_dif[i]))
    std.append(np.nanstd(abs_dif[i]))
    mean.append(np.nanmean(abs(abs_dif[i])))