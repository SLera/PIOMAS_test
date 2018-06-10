#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 13:56:10 2018

@author: valeria
"""
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
#import datetime
#from scipy import stats
import sys
sys.path.append('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')
import Vi_functions as V



INDIR = '/home/valeria/NIERSC/Projects/Bashmachnikov/data/PIOMAS/'
OUTDIR = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test/Comparison_obs/'
months_winter = np.array([10,11,12,1,2,3,4])
months_dc = np.array([1,2,3,4])

##PIOMAS
fname_PIOMAS = INDIR+'PIOMAS_vol_export.txt'
outname_PIOMAS = OUTDIR+'txt/'+'Monthly_mean_ViFlix_std_PIOMAS_1979-2016_OctApr.txt'
Vif_PIOMAS, dates_PIOMAS = V.read_PIOMAS_Vi_flux(fname_PIOMAS)
years_PIOMAS = np.arange(1978,2017)
##winter
monthly_Vif_PIOMAS, mean_monthly_Vif_PIOMAS, std_monthly_Vif_PIOMAS = V.calculate_monthly_stats_Vif(years_PIOMAS, months_winter, Vif_PIOMAS, dates_PIOMAS, outname_PIOMAS)
###dc
#outname_PIOMAS = 'Monthly_mean_ViFlix_std_PIOMAS_1979-2016_JanApr.txt'
#monthly_Vif_PIOMAS_dc, mean_monthly_Vif_PIOMAS_dc, std_monthly_Vif_PIOMAS_dc = calculate_monthly_stats_Vif(years_PIOMAS, months_dc, Vif_PIOMAS, dates_PIOMAS, outname_PIOMAS)

##Cryosat
fname_Cryosat = INDIR+'CRYOSAT_vol_export.txt'
outname_Cr = OUTDIR+'txt/'+'Monthly_mean_ViFlix_std_Cryosat_2010-2017_OctApr.txt'
Vif_Cr, dates_Cr = V.read_Cryosat_Vi_flux(fname_Cryosat)
years_Cr = np.arange(2010,2018,1)
##winter
monthly_Vif_Cr, mean_monthly_Vif_Cr, std_monthly_Vif_Cr = V.calculate_monthly_stats_Vif(years_Cr, months_winter, Vif_Cr, dates_Cr, outname_Cr)
###dc
#outname_Cr = 'Monthly_mean_ViFlix_std_Cryosat_2010-2017_JanApr.txt'
#Vif_period_Cr_dc, dates_period_Cr_dc, mean_monthly_Vif_Cr_dc, std_monthly_Vif_Cr_dc = calculate_monthly_stats_Vif(years_Cr, months_dc, Vif_Cr, dates_Cr, outname_Cr)

#Kwok_2004
fname_k04=INDIR+'Kwok_vol_export_2004.txt'
outname_k04 =OUTDIR+'txt/'+ 'Monthly_mean_ViFlix_std_Kwok_1991-1999_OctApr.txt'
Vif_k04, dates_k04 = V.read_Kwk_Sprn_Vi_flux(fname_k04)
years_k04 = np.arange(1991,2000,1)
##winter
monthly_Vif_k04, mean_monthly_Vif_k04, std_monthly_Vif_k04 = V.calculate_monthly_stats_Vif(years_k04, months_winter, Vif_k04, dates_k04, outname_k04)

##    
##import Spreen_2009 volume flux
fname_s09=INDIR+'Spreen_vol_export_2009.txt'
outname_s09 = OUTDIR+'txt/'+'Monthly_mean_ViFlix_std_Spreen_2003-2008_OctApr.txt'
Vif_s09, dates_s09 = V.read_Kwk_Sprn_Vi_flux(fname_s09)
years_s09 = np.arange(2003,2009,1)
##winter
monthly_Vif_s09, mean_monthly_Vif_s09, std_monthly_Vif_s09 = V.calculate_monthly_stats_Vif(years_s09, months_winter, Vif_s09, dates_s09, outname_s09)


#PLOTTING
plt.figure()
#PIOMAS
plt.plot(mean_monthly_Vif_PIOMAS,  c = 'royalblue', label = 'PIOMAS')
yerr=std_monthly_Vif_PIOMAS
plt.errorbar(range(len(yerr)),mean_monthly_Vif_PIOMAS,yerr,lw=0, ecolor='royalblue', elinewidth=1, capsize=5)

plt.plot(mean_monthly_Vif_Cr, c = 'darkorange', label = 'Cryosat2/SMOS')
yerr=std_monthly_Vif_Cr
plt.errorbar(range(len(yerr)),mean_monthly_Vif_Cr,yerr,lw=0, ecolor='darkorange', elinewidth=1, capsize=5)

plt.plot(mean_monthly_Vif_k04, c = 'plum',  label = 'Kwok2004')
yerr=std_monthly_Vif_k04
plt.errorbar(range(len(yerr)),mean_monthly_Vif_k04,yerr,lw=0, ecolor='plum', elinewidth=1, capsize=5)

plt.plot(mean_monthly_Vif_s09, c = 'mediumspringgreen',  label = 'Spreen2009')
yerr=std_monthly_Vif_s09
plt.errorbar(range(len(yerr)),mean_monthly_Vif_s09,yerr,lw=0, ecolor='mediumspringgreen', elinewidth=1, capsize=5)

plt.xticks(np.arange(len(months_winter)),[str(m) for m in months_winter])
plt.xlabel('Months')
plt.ylabel('km^3')
plt.title('Mean monthly sea ice volume flux')
plt.grid()
plt.tight_layout()
plt.legend()
plt.show()
plt.savefig(OUTDIR+'figs/'+'PIOMAS_Viflux_comparison_cycle.pdf')
#
