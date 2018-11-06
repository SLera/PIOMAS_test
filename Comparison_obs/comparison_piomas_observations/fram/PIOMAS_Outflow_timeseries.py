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
import sys
sys.path.append('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')
import Vi_functions as V

##PLOT SETTINGS
##Create colormap
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

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

##I/O 
INDIR = '/home/lera/NIERSC/Projects/Bashmachnikov/data/PIOMAS/'
OUTDIR = '/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/POLAR2018/'
# import PIOMAS volume flux
fname_PIOMAS = INDIR + 'piomas_vol_export_v2.txt'
years_un = np.arange(1979, 2017, 1)
months_un = np.arange(1, 13, 1)
Vif_PIOMAS, dates_PIOMAS = V.read_PIOMAS_Vi_flux(fname_PIOMAS)

# import Cryosat volume flux
years_cr = np.arange(2010, 2017, 1)
fname_Cryosat = INDIR + 'CRYOSAT_vol_export.txt'
Vif_Cr, dates_Cr = V.read_Cryosat_Vi_flux(fname_Cryosat)

# import Kwok_2004 volume flux
fname_k04 = INDIR + 'Kwok_vol_export_2004.txt'
Vif_k04, dates_k04 = V.read_Kwk_Sprn_Vi_flux(fname_k04)

# import Spreen_2009 volume flux
fname_s09 = INDIR + 'Spreen_vol_export_2009.txt'
Vif_s09, dates_s09 = V.read_Kwk_Sprn_Vi_flux(fname_s09)

# plot timeseries Riker, PIOMAS, Kwok04, Spreen09
years = np.array([int(d.year) for d in dates_PIOMAS])
#start PIOMAS time series from year 1991
ind = np.where(np.array(years) == 1991)
i = ind[0][0]

#set fig size in inch ((x,y))
fsize = cm2inch((19,8.2))
plt.figure(figsize=fsize)

plt.plot(dates_PIOMAS[i:], Vif_PIOMAS[i:],  lw = 1.7, c="#3F5D7D", label='PIOMAS')
plt.plot(dates_k04, Vif_k04, c=tableau20[12], lw = 1.5, alpha = 0.7,label='Kwok et al. (2004)')
plt.plot(dates_s09, Vif_s09, c=tableau20[4], lw = 1.5,alpha = 0.7, label='Spreen et al. (2009)')
plt.plot(dates_Cr[:-4], Vif_Cr[:-4], c=tableau20[2], lw = 1.5,alpha = 0.7,  label='Ricker et al. (2017)')


plt.title('Monthly ice volume exprort through Fram Strait', fontsize = 8)
plt.ylabel('km^3', fontsize = 8)
plt.ylim((-60,650))
plt.xlim((dates_PIOMAS[i], datetime.datetime(2017,4,1)))
plt.xticks( fontsize = 8)
plt.yticks( fontsize = 8)
plt.legend(ncol = 2, fontsize = 6 )
plt.grid(axis = 'y', linestyle='--')
plt.tight_layout()
# plt.show()
plt.savefig(OUTDIR+'PIOMAS_Viflux_comparison_ts.pdf')
plt.savefig(OUTDIR+'PIOMAS_Viflux_comparison_ts.png')
plt.close()
