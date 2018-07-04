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
INDIR = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/'
OUTDIR = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/figs/Greenland/Hi_mean_monthly/'

chi_months = np.load(INDIR+'CryosatHi_Greenland_sorted_monthly')
chi_u_months = np.load(INDIR+'CryosatUnc_Greenland_sorted_monthly')
phi_months = np.load(INDIR+'PIOMAS_Greenland_sorted_monthly')

clat = np.load(INDIR+'clat_Greenland')
clon = np.load(INDIR+'clon_Greenland')
months = np.array([10,11,12,1,2,3,4])

for m in range(len(months)):
    print m
    c = chi_months[m]
    p = phi_months[m]
    cu = chi_u_months[m]
    e,f,g = np.shape(chi_months[m])
    mean_c = np.zeros(e)
    mean_p = np.zeros(e)
    monthly_dif = np.zeros((f,g))
    for i in range(f):
        for j in range(g):
            for h in range(e):
                value_c = c[h][i,j]
                mean_c[h] = value_c
                value_p = p[h][i,j]
                mean_p[h] = value_p
            monthly_dif[i,j]=np.nanmean(mean_p-mean_c)
    ind = np.where(monthly_dif==0)
    monthly_dif[ind]=np.nan
    outfname = OUTDIR+'Monthly_mean_dif_'+str(months[m])+'.png'
    title = 'PIOMAS-Cryosat, m'
    Read_Hi_data.plot_map_Greenland(monthly_dif, (-5,5), clon, clat, outfname, title)