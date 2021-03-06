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
from mpl_toolkits.basemap import Basemap
import datetime
from scipy import stats
import sys
sys.path.append('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')
import Read_Hi_data


def extract_date_PIOMAS(path):
    return datetime.datetime.strptime(path[-9:-3], '%Y%m')


def extract_date_Cr(path):
    return datetime.datetime.strptime(path[-15:-9], '%Y%m')


def find_absdif(Phi, Chi):
    param = Phi - Chi
    return param

def filter_nan(array1,array2):
    array1_f1 = array1[~np.isnan(array1)]
    array2_f1 = array2[~np.isnan(array1)]
    array1_f2 = array1_f1[~np.isnan(array2_f1)]
    array2_f2 = array2_f1[~np.isnan(array2_f1)]
    return array1_f2, array2_f2

# GREATE COLORMAP
# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

OUTDIR = '/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/figs/Greenland/'

# IMPORT CORRESPONDING ARRAYS OF Hi FROM PIOMAS AND CRYOSAT (Create_PIOMASvsCryosat_Hiarrays.py)

months = np.array([10,11,12,1,2,3,4])

chi_months = np.load('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/CryosatHi_full_sorted_monthly')
phi_months = np.load('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/PIOMAS_full_sorted_monthly')
chi_u_months = np.load('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/CryosatUnc_full_sorted_monthly')
clat = np.load('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/clat_full')
clon = np.load('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/clon_full')

# PLOT MONTHLY SCATTERPLOTS
fig, axes = plt.subplots(figsize=(7.38,4), nrows=2, ncols=4, sharex=True, sharey=True)
fig.subplots_adjust(wspace=0.05, hspace=0.20, top=0.95, bottom=0.05)

p_mean = []
c_mean = []
c_u_mean = []

for m in range(len(months)):
    print m
    cr_m = []
    p_m = []
    plt.subplot('24' + str(m + 1))
    plt.title(str(months[m]))
    for i in range(len(chi_months[m])):
        c = chi_months[m][i]
        c_u = chi_u_months[m][i]
        p = phi_months[m][i]
        cr, ps = filter_nan(c, p)
        cr_u, _ = filter_nan(c_u, p)
        # exclude Cryosat Hi with uncertainties > measurements
        ps = ps[np.where(cr>cr_u)]
        cr = cr[np.where(cr>cr_u)]
        p_mean.append(ps.mean())
        c_mean.append(cr.mean())
        c_u_mean.append(cr_u.mean())
        plt.scatter(cr,ps, s = 5, facecolor='none',alpha = 0.2,edgecolor = tableau20[m],label = str(m))
        # plt.errorbar(cr,p, xerr=cr_u/2, color = tableau20[m], linestyle = 'None')
        cr_m.extend(cr)
        p_m.extend(ps)

        test = np.where(ps<0)
        if len(test[0]!=0):
            print "negative Hi,", m, i

    s, i, r, p, std = stats.linregress(np.array(cr_m), np.array(p_m))
    plt.text(0.5,7.1, 'slope = '+ str(np.round(s,2)))
    plt.text(0.5,6.3, 'r = '+ str(np.round(r,2)))
    plt.xlim((0,8))
    plt.ylim((0,8))
    plt.plot([-1,9],[-1,9],'k-', lw = 0.5)
    plt.xticks(np.arange(0,9,1))
    plt.yticks(np.arange(0,9,1))
    # plt.axes().set_aspect('equal')

plt.subplot(248)
plt.title('Mean monthly')
plt.scatter(c_mean,p_mean,s = 5,facecolor='black',edgecolor = 'black')
# plt.errorbar(cr,p, xerr=cr_u/2, color = tableau20[m], linestyle = 'None')
slope2, intercept2, r_value2, p_value2, std_er2 = stats.linregress(c_mean,p_mean)
plt.text(1.2, 2.83, 'slope = ' + str(np.round(slope2, 2)))
plt.text(1.2, 2.6, 'r = ' + str(np.round(r_value2, 2)))
plt.xlim((1, 3))
plt.ylim((1, 3))
plt.plot([-1, 3], [-1, 3], 'k-', lw=0.5)
plt.xticks(np.arange(1, 3.1, 0.5))
plt.yticks(np.arange(1, 3.1, 0.5))

plt.tight_layout()
plt.show()
plt.savefig(OUTDIR+'Monthly_Hi_scatter_PIOMASvsCryosat2.png')
plt.close()
