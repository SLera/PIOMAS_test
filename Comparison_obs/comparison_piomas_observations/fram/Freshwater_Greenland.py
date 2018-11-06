#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 13:56:10 2018

@author: valeria
"""
import numpy as np
import matplotlib
matplotlib.use('qt5agg')
from matplotlib import pyplot as plt
import datetime
from scipy import stats
import glob, os
import sys
sys.path.append('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')
import Vi_functions as V
import Read_Hi_data


fname_Gr = '/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/output/Greenland/EASE2N_grid_v2/PIOMAS_Vi_Greenland_EASE2N.txt'
years,months,vi = np.loadtxt(fname_Gr, unpack=True)
dates = [datetime.date(np.int(years[i]), np.int(months[i]),1) for i in range(len(years))]


# import PIOMAS volume flux
INDIR_outflow = '/home/lera/NIERSC/Projects/Bashmachnikov/data/PIOMAS/'
fname_Fram = INDIR_outflow + 'piomas_vol_export_v2.txt'
fname_66 = INDIR_outflow + 'piomas_vol_export_66.txt'
Vif_Fram, dates_Fram = V.read_PIOMAS_Vi_flux(fname_Fram)
Vif_66, dates_66 = V.read_PIOMAS_Vi_flux(fname_66)
# OUTDIR = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/figs/Greenland/'
V_melt = []

for i in range(len(vi)-1):
    print i
    v_m = vi[i]
    v_m1 = vi[i+1]
    fram = Vif_Fram[i]
    bb = Vif_66[i]
    melt = v_m1-v_m+fram-bb
    print melt
    V_melt.append(melt)

print np.shape(np.array(years[1:])), np.shape(np.array(months[1:])), np.shape(np.array(V_melt))


table = np.column_stack((np.array(years[1:]),np.array(months[1:]),np.array(V_melt)))
np.savetxt(INDIR_outflow+'PIOMAS_Vi_freeze_melt.txt', table,  header = 'timestamp, Sea ice volume freeze(+) or melt (-) in km^3', fmt='%1.2f')
