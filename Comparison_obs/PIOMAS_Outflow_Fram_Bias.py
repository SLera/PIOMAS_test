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

def bias_Vi(model,obs):
    n = len(model)
    b = obs-model
#    for i in range(len(dataset1)):
#        b[i]=dataset1[i]-dataset2[i]read_PIOMAS_Vi_flux()
    bias = np.sum(b)/n
    return bias

def rmse_Vi(model,obs):
    n = len(model)
    r = np.zeros((n))
    for i in range(len(model)):
        r[i]=np.sqrt((obs[i]-model[i])**2)
    rmse = np.sum(r/n)
    return rmse

def rpd_Vi(model,obs):
    n = len(model)
    dif = obs-model
    rdif = dif/model
    rdif_sum = np.sum(rdif)
    rpd = rdif_sum/n
    return rpd*100


INDIR = '/home/lera/NIERSC/Projects/Bashmachnikov/data/PIOMAS/'
OUTDIR = '/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/'
months_winter = np.array([10,11,12,1,2,3,4])
months_dc = np.array([1,2,3,4])

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

# find corresponding data for POIMAS/Riker
dates_pm_cr = []
Viflux_pm_cr = []
winter_months = [10, 11, 12, 1, 2, 3, 4]

for n in range(0, len(dates_Cr[:-4])):
    date_cr = dates_Cr[n]
    if date_cr.month in winter_months:
        timelag = dates_Cr[n] - np.array(dates_PIOMAS)
        timelag_ind = np.where(abs(timelag) == abs(timelag).min())[0][0]
        test = timelag[timelag_ind]
        if test.days != 0:
            print test.days, 'date_cr:', date_cr, 'date_pm:', dates_PIOMAS[timelag_ind]
        else:
            dates_pm_cr.append(dates_PIOMAS[timelag_ind])
            Viflux_pm_cr.append(Vif_PIOMAS[timelag_ind])

mask = ~np.isnan(Vif_Cr)

bias_Cr = bias_Vi(Viflux_pm_cr,Vif_Cr[mask][:-4])
rmse_Cr = rmse_Vi(Viflux_pm_cr,Vif_Cr[mask][:-4])
rpd_Cr = rpd_Vi(Viflux_pm_cr,Vif_Cr[mask][:-4])
print 'bias_cr', bias_Cr
print 'rmse_cr', rmse_Cr
print 'rpd_cr', rpd_Cr

# find corresponding data for scatterplot POIMAS/Kwok2004
dates_pm_k04_winter = []
Viflux_pm_k04_winter = []
dates_k04_winter = []
Viflux_k04_winter = []

winter_months = [10, 11, 12, 1, 2, 3, 4, 5]

for n in range(0, len(dates_k04)):
    date_k04 = dates_k04[n]
    if date_k04.month in winter_months:
        dates_k04_winter.append(date_k04)
        Viflux_k04_winter.append(Vif_k04[n])
        timelag = dates_k04[n] - np.array(dates_PIOMAS)
        timelag_ind = np.where(abs(timelag) == abs(timelag).min())[0][0]
        test = timelag[timelag_ind]
        if test.days != 0:
            print test.days, 'date_k04:', date_k04, 'date_pm:', dates_PIOMAS[timelag_ind]
        else:
            dates_pm_k04_winter.append(dates_PIOMAS[timelag_ind])
            Viflux_pm_k04_winter.append(Vif_PIOMAS[timelag_ind])

bias_k04 = bias_Vi(np.array(Viflux_pm_k04_winter), np.array(Viflux_k04_winter))
rmse_k04 = rmse_Vi(np.array(Viflux_pm_k04_winter), np.array(Viflux_k04_winter))
rpd_k04 = rpd_Vi(np.array(Viflux_pm_k04_winter), np.array(Viflux_k04_winter))
print 'bias_k04', bias_k04
print 'rmse_k04', rmse_k04
print 'rpd_k04', rpd_k04

# create data for scatterplot POIMAS/Spreen09
dates_pm_s09 = []
Viflux_pm_s09 = []
winter_months = [10, 11, 12, 1, 2, 3, 4]

for n in range(0, len(dates_s09)):
    date_s09 = dates_s09[n]
    if date_s09.month in winter_months:
        timelag = dates_s09[n] - np.array(dates_PIOMAS)
        timelag_ind = np.where(abs(timelag) == abs(timelag).min())[0][0]
        test = timelag[timelag_ind]
        if test.days != 0:
            print test.days, 'date_s09:', date_s09, 'date_pm:', dates_PIOMAS[timelag_ind]
        else:
            dates_pm_s09.append(dates_PIOMAS[timelag_ind])
            Viflux_pm_s09.append(Vif_PIOMAS[timelag_ind])

mask_s09 = ~np.isnan(Vif_s09)
slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(Vif_s09[mask_s09], Viflux_pm_s09)

bias_s09 = bias_Vi(np.array(Viflux_pm_s09), np.array(Vif_s09[mask_s09]))
rmse_s09 = rmse_Vi(np.array(Viflux_pm_s09), np.array(Vif_s09[mask_s09]))
rpd_s09= rpd_Vi(np.array(Viflux_pm_s09), np.array(Vif_s09[mask_s09]))
print 'bias_s09', bias_s09
print 'rmse_s09', rmse_s09
print 'rpd_s09', rpd_s09
