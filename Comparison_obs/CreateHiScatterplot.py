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

chi_months = [[] for x in xrange(len(months))]
chi_u_months = [[] for x in xrange(len(months))]
phi_months = [[] for x in xrange(len(months))]

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

                # chi = Read_Hi_data.extract_Greenland(chi)
                # chi_u = Read_Hi_data.extract_Greenland(chi_u)
                # phi = Read_Hi_data.extract_Greenland(phi)

                # phi = Read_Hi_data.cut_domain(phi)
                # clat = Read_Hi_data.cut_domain(clat)
                # clon = Read_Hi_data.cut_domain(clon)
                # chi_u = Read_Hi_data.cut_domain(chi_u)
                # chi = Read_Hi_data.cut_domain(chi)
                ind_month = np.where(months == date_PIOMAS.month)[0][0]
                chi_months[ind_month].append(chi)
                chi_u_months[ind_month].append(chi_u)
                phi_months[ind_month].append(phi)

c = np.array(chi_months)
c.dump('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/CryosatHi_full_sorted_monthly')
p = np.array(phi_months)
p.dump('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/PIOMAS_full_sorted_monthly')
cu = np.array(chi_u_months)
cu.dump('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/CryosatUnc_full_sorted_monthly')
clat.dump('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/clat_full')
clon.dump('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/clon_full')
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

# PLOT MONTHLY SCATTERPLOTS
plt.figure()
for m in range(len(months)):
    print m
    plt.figure()
    for i in range(len(chi_months[m])):
        c = chi_months[m][i]
        cr=c[~np.isnan(c)]
        cr_u = chi_u_months[m][i]
        cr_u=cr_u[~np.isnan(c)]
        cr = cr - cr_u
        p = phi_months[m][i]
        p=p[~np.isnan(c)]
        print np.shape(cr), np.shape(p)
        plt.scatter(cr,p,facecolor='none',alpha = 0.05,edgecolor = tableau20[m],label = str(m))
        # plt.errorbar(cr,p, xerr=cr_u/2, color = tableau20[m], linestyle = 'None')
        plt.xlim((0,10))
        plt.ylim((0,10))
plt.legend()
plt.show()

# PLOT MEAN MONTHLY SCATTER
plt.figure()
p_mean = []
c_mean = []
c_u_mean = []
for m in range(len(months)):
    print m
    for i in range(len(chi_months[m])):
        c = chi_months[m][i]
        cr=c[~np.isnan(c)]
        cr_u = chi_u_months[m][i]
        cr_u=cr_u[~np.isnan(c)]
        p = phi_months[m][i]
        p=p[~np.isnan(c)]
        p_mean.append(p.mean())
        c_mean.append(cr.mean())
        c_u_mean.append(cr_u.mean())

c_mean = np.array(c_mean)
c_u_mean = np.array(c_u_mean)
c_mean = c_mean - c_u_mean
p_mean = np.array(p_mean)
c_mean = c_mean[~np.isnan(p_mean)]
p_mean = p_mean[~np.isnan(p_mean)]

plt.scatter(c_mean,p_mean,facecolor='none',edgecolor = tableau20[m])
# plt.errorbar(cr,p, xerr=cr_u/2, color = tableau20[m], linestyle = 'None')
slope2, intercept2, r_value2, p_value2, std_er2 = stats.linregress(c_mean,p_mean)
line2 = slope2*np.array(c_mean)+intercept2

print 'PIOMASvsCryosat_mean: r, slope, p',  r_value2, slope2, p_value2
print 'variance PIOMAS_mean:', np.var(p_mean)
print 'variance Cryosat_mean:', np.var(c_mean)
plt.plot(c_mean,line2,ls ='-', c = tableau20[m],lw=0.9)
plt.plot()
plt.xlim((0.5,3))
plt.ylim((0.5,3))
plt.legend()
plt.show()

