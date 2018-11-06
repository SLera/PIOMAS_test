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
import sys
from scipy import stats
sys.path.append('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')
import Read_Hi_data


INDIR = '/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/'
OUTDIR = '/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/figs/Greenland/'

chi_months = np.load(INDIR+'CryosatHi_full_sorted_monthly')
chi_u_months = np.load(INDIR+'CryosatUnc_full_sorted_monthly')
phi_months = np.load(INDIR+'PIOMAS_full_sorted_monthly')
clat = np.load(INDIR+'clat_full')
clon = np.load(INDIR+'clon_full')

months = np.array([10,11,12,1,2,3,4])

#DOWNSCALE
# chi_months_sc = [[] for x in xrange(len(months))]
#
# phi_months_sc = [[] for x in xrange(len(months))]
#
#
# for m in range (len(months)):
#     scale = 20
#     for i in range(len(chi_months[m])):
#         chi_months_sc[m].append(Read_Hi_data.downscale_array_running_mean(chi_months[m][i],scale))
#         phi_months_sc[m].append(Read_Hi_data.downscale_array_running_mean(phi_months[m][i], scale))
#
# chi_months = chi_months_sc
# phi_months = phi_months_sc
#
# all_c_months = []
# all_p_months = []

# # MONTHLY CORRELATION
# plt.figure()
# for m in range(len(months)):
#     print m
#     plt.figure()
#     for i in range(len(chi_months[m])):
#         c = chi_months[m][i]
#         cr=c[~np.isnan(c)]
#         p = phi_months[m][i]
#         p=p[~np.isnan(c)]
#         plt.scatter(cr,p,facecolor='none',alpha = 0.5,edgecolor = tableau20[m],label = str(m))
#         all_c_months.extend(cr.flatten())
#         all_p_months.extend(p.flatten())
#         plt.xlim((0,10))
#         plt.ylim((0,10))
#     plt.legend()
#     r,p = stats.pearsonr(np.array(all_c_months),np.array(all_p_months))
#     print months[m], r,p
# plt.show()

#GRID CELL CORRELATION
min_sample_len = 10 #min number of observations in a grid cell for correlation calc
cor_array = np.zeros(np.shape(chi_months[0][0]))
cor_array.fill(np.nan)
for i in range(np.shape(cor_array)[0]):
    print i, 'from', np.shape(cor_array)[0]
    for j in range(np.shape(cor_array)[1]):
        values_c = []
        values_p = []
        for m in range(len(months)):
            c = chi_months[m]
            p = phi_months[m]
            for l in range(len(c)):
                values_c.append(c[l][i,j])
                values_p.append(p[l][i,j])
            if not all( np.isnan(v) for v in values_c):
                values_c_n = np.array(values_c)[~np.isnan(values_c)]
                values_p_n = np.array(values_p)[~np.isnan(values_c)]
                values_c_n = values_c_n[values_p_n>=0]
                values_p_n = values_p_n[values_p_n>=0]
                if len(values_c_n)>min_sample_len:
                    r,p = stats.pearsonr(values_c_n,values_p_n)
                    #if p < 0.05:
                    cor_array[i,j]= r

    # ind = np.where(monthly_dif==0)
    # monthly_dif[ind]=np.nan
outfname = OUTDIR+'Correlation_map_Arctic'+'.pdf'
title = 'PIOMASvsCryosat,cor_coef'
Read_Hi_data.plot_map_Greenland(cor_array, (0,0.8), clon, clat, outfname, title)

outfname = OUTDIR+'Correlation_map_Greenland'+'.pdf'
title = 'PIOMASvsCryosat,cor_coef'
Read_Hi_data.plot_map_Greenland(Read_Hi_data.cut_domain(cor_array), (0,0.8), Read_Hi_data.cut_domain(clon), Read_Hi_data.cut_domain(clat), outfname, title)

c = Read_Hi_data.extract_Greenland(cor_array)
print 'Greenland Sea mean cor corf', np.nanmean(c)
print 'Greenland Sea median cor corf', np.nanmedian(c)
