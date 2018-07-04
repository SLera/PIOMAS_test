#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 13:56:10 2018

@author: valeria
"""
import numpy as np
from matplotlib import pyplot as plt
import datetime
from scipy import stats
import sys

sys.path.append('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')
import Vi_functions as V

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

INDIR = '/home/valeria/NIERSC/Projects/Bashmachnikov/data/PIOMAS/'
OUTDIR = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/POLAR2018/'
# import PIOMAS volume flux
fname_PIOMAS = INDIR + 'PIOMAS_vol_export.txt'
years_un = np.arange(1979, 2017, 1)
months_un = np.arange(1, 13, 1)
Vif_PIOMAS, dates_PIOMAS = V.read_PIOMAS_Vi_flux(fname_PIOMAS)

# c = 1.710949376507358
c = 1
Vif_PIOMAS, dates_PIOMAS = V.read_PIOMAS_Vi_flux(fname_PIOMAS)


ys = np.arange(1979, 2017)
mean_season_p = np.zeros(len(ys))
months = [1,2,3,4,10,11,12]
for y in range(len(ys)):
    list_flux_p = []
    for i in range(len(dates_PIOMAS)):
        date = dates_PIOMAS[i]
        if ys[y]==date.year and date.month in months:
            list_flux_p.append(Vif_PIOMAS[i])
    mean_season_p[y] = np.nanmean(np.array(list_flux_p))

s_p, i_p, r_p, p_p, std_p =  stats.linregress(ys, mean_season_p)

print s_p, i_p, r_p, p_p, std_p