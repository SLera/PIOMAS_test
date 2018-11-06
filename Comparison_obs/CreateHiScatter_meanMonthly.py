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

# PLOT MEAN MONTHLY SCATTER
plt.figure()
p_mean = []
c_mean = []
c_u_mean = []
p_mean_filt = []
c_mean_filt = []

for m in range(len(months)):
    print m
    for i in range(len(chi_months[m])):
        c = chi_months[m][i]
        c_u = chi_u_months[m][i]
        p = phi_months[m][i]
        cr, ps = filter_nan(c, p)
        cr_u, _ = filter_nan(c_u, p)
        # exclude Cryosat Hi with uncertainties > measurements
        ps_filt = ps[np.where(cr > cr_u)]
        cr_filt = cr[np.where(cr > cr_u)]
        p_mean.append(ps.mean())
        c_mean.append(cr.mean())
        c_u_mean.append(cr_u.mean())
        p_mean_filt.append(ps_filt.mean())
        c_mean_filt.append(cr_filt.mean())

c_mean = np.array(c_mean)
p_mean = np.array(p_mean)
# c_mean = c_mean[~np.isnan(p_mean)]
# p_mean = p_mean[~np.isnan(p_mean)]

c_mean_filt = np.array(c_mean_filt)
p_mean_filt = np.array(p_mean_filt)
# c_mean_filt = c_mean_filt[~np.isnan(p_mean_filt)]
# p_mean_filt = p_mean_filt[~np.isnan(p_mean_filt)]


plt.scatter(c_mean,p_mean,facecolor='none',edgecolor = 'b', label = 'all Hi_CS')
plt.scatter(c_mean_filt, p_mean_filt, facecolor='none', edgecolor='r', label='Hi_CS>Hi_CS_uncertainty')
# plt.errorbar(cr,p, xerr=cr_u/2, color = tableau20[m], linestyle = 'None')
slope2, intercept2, r_value2, p_value2, std_er2 = stats.linregress(c_mean,p_mean)
line2 = slope2*np.array(c_mean)+intercept2

slope3, intercept3, r_value3, p_value3, std_er3 = stats.linregress(c_mean_filt, p_mean_filt)
line2 = slope2 * np.array(c_mean) + intercept2
print 'PIOMASvsCryosat_mean: r, slope, p',  r_value2, slope2, p_value2
print 'variance PIOMAS_mean:', np.var(p_mean)
print 'variance Cryosat_mean:', np.var(c_mean)

print 'PIOMASvsCryosat_mean -UN: r, slope, p', r_value3, slope3, p_value3
print 'variance PIOMAS_mean -UN:', np.var(p_mean_filt)
print 'variance Cryosat_mean -UN:', np.var(c_mean_filt)
# plt.plot(c_mean,line2,ls ='-', c = tableau20[m],lw=0.9)
plt.plot()
plt.xlabel('CS, m')
plt.ylabel('PIOMAS, m')
plt.xlim((0.5,3))
plt.ylim((0.5,3))
plt.show()

