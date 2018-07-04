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

Vif_PIOMAS_c, dates_PIOMAS = V.read_PIOMAS_Vi_flux(fname_PIOMAS)

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

# plot timeseries Cryosat, PIOMAS, Kwock04, Spreen09
years = np.array([int(d.year) for d in dates_PIOMAS])
ind = np.where(np.array(years) == 1991)
i = ind[0][0]

from matplotlib import font_manager
font_manager._rebuild()
import matplotlib
matplotlib.rcParams['font.family'] = 'Montserrat'

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

fsize = cm2inch((23,10))
plt.figure(figsize=fsize)

plt.plot(dates_PIOMAS[i:], Vif_PIOMAS[i:],  lw = 1.7, c="#3F5D7D", label='PIOMAS')
# plt.plot(dates_PIOMAS[i:], Vif_PIOMAS_c[i:], ls='--', c=tableau20[0], label='PIOMAS')
plt.plot(dates_k04, Vif_k04, c=tableau20[12], lw = 1.5, alpha = 0.7,label='Kwok et al. (2004)')
plt.plot(dates_s09, Vif_s09, c=tableau20[4], lw = 1.5,alpha = 0.7, label='Spreen et al. (2009)')
plt.plot(dates_Cr[:-4], Vif_Cr[:-4], c=tableau20[2], lw = 1.5,alpha = 0.7,  label='Ricker et al. (2017)')

plt.title('Monthly ice volume export though Fram Strait', fontsize = 20)
plt.ylabel('km^3', fontsize = 18)
plt.xlabel('Observations', fontsize = '18')
plt.ylim((-60,650))
plt.xticks( fontsize = 16)
plt.yticks( fontsize = 16)
plt.legend(ncol = 2, fontsize = 12 )
plt.grid(axis = 'y', linestyle='--')
plt.tight_layout()
# plt.show()
plt.savefig(OUTDIR+'PIOMAS_Viflux_comparison_ts.pdf')
plt.close()
# create data for scatterplot POIMAS/CRYOSAT
Vif_PIOMAS = Vif_PIOMAS_c

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
slope, intercept, r_value, p_value, std_err = stats.linregress(Vif_Cr[mask][:-4], Viflux_pm_cr)
line = slope * Vif_Cr[mask][:-4] + intercept
print 'PIOMASvsCryosat:', r_value, slope, p_value
print 'variance PIOMAS:', np.var(Viflux_pm_cr)
print 'variance Cryosat:', np.var(Vif_Cr[mask][:-4])


ys = np.arange(2010, 2017)
mean_season = np.zeros(len(ys))
mean_season_p = np.zeros(len(ys))

for y in range(len(ys)):
    list_flux_p = []
    list_flux = []
    for i in range(len(dates_pm_cr)):
        date = dates_pm_cr[i]
        if ys[y]==date.year:
            list_flux_p.append(Viflux_pm_cr[i])
            list_flux.append(Vif_Cr[i])
    mean_season[y] = np.nanmean(np.array(list_flux))
    mean_season_p[y] = np.nanmean(np.array(list_flux_p))

s, i, r, p, std =  stats.linregress(ys, mean_season)
s_p, i_p, r_p, p_p, std_p =  stats.linregress(ys, mean_season_p)

print s, r**2
print s_p, r_p**2

# create data for scatterplot POIMAS/Kwok2004
dates_pm_k04_winter = []
Viflux_pm_k04_winter = []
dates_k04_winter = []
Viflux_k04_winter = []

dates_pm_k04_summer = []
Viflux_pm_k04_summer = []
dates_k04_summer = []
Viflux_k04_summer = []

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
            print test.days, 'date_k04:', date_cr, 'date_pm:', dates_PIOMAS[timelag_ind]
        else:
            dates_pm_k04_winter.append(dates_PIOMAS[timelag_ind])
            Viflux_pm_k04_winter.append(Vif_PIOMAS[timelag_ind])
    else:
        dates_k04_summer.append(date_k04)
        Viflux_k04_summer.append(Vif_k04[n])
        timelag = dates_k04[n] - np.array(dates_PIOMAS)
        timelag_ind = np.where(abs(timelag) == abs(timelag).min())[0][0]
        test = timelag[timelag_ind]
        if test.days != 0:
            print test.days, 'date_k04:', date_cr, 'date_pm:', dates_PIOMAS[timelag_ind]
        else:
            dates_pm_k04_summer.append(dates_PIOMAS[timelag_ind])
            Viflux_pm_k04_summer.append(Vif_PIOMAS[timelag_ind])

slope1, intercept1, r_value1, p_value1, std_er1 = stats.linregress(Viflux_k04_winter, Viflux_pm_k04_winter)
line1 = slope1 * np.array(Viflux_k04_winter) + intercept1
print 'PIOMASvsKwok04:', r_value1, slope1, p_value1
print 'variance PIOMAS:', np.var(Viflux_pm_k04_winter)
print 'variance Kwok04:', np.var(Viflux_k04_winter)

# ys = np.arange(1991, 2000)
# mean_season = np.zeros(len(ys))
# mean_season_p = np.zeros(len(ys))
# for y in range(len(ys)):
#     list_flux_p = []
#     list_flux = []
#     for i in range(len(dates_pm_k04_winter)):
#         date = dates_pm_k04_winter[i]
#         if ys[y]==date.year:
#             list_flux_p.append(Viflux_pm_k04_winter[i])
#             list_flux.append(Viflux_k04_winter[i])
#     mean_season[y] = np.nanmean(np.array(list_flux))
#     mean_season_p[y] = np.nanmean(np.array(list_flux_p))
#
# s, i, r, p, std =  stats.linregress(ys, mean_season)
# s_p, i_p, r_p, p_p, std_p =  stats.linregress(ys, mean_season_p)
#
# print s, r**2
# print s_p, r_p**2

##find summer_mean
years_k04_summer = np.arange(1992, 2000)
Viflux_k04_summer_mean = np.zeros(8)
Viflux_pm_k04_summer_mean = np.zeros(8)

for i in range(len(Viflux_k04_summer)):
    date_k04 = dates_k04_summer[i]
    for j in range(len(years_k04_summer)):
        if date_k04.year == years_k04_summer[j]:
            Viflux_k04_summer_mean[j] += Viflux_k04_summer[i]
            Viflux_pm_k04_summer_mean[j] += Viflux_pm_k04_summer[i]

slope2, intercept2, r_value2, p_value2, std_er2 = stats.linregress(Viflux_k04_summer_mean, Viflux_pm_k04_summer_mean)
line2 = slope2 * np.array(Viflux_k04_summer_mean) + intercept2
print 'PIOMASvsKwok04_Jun-Sep_mean:', r_value2, slope2, p_value2
print 'variance PIOMAS_Jun-Sep_mean:', np.var(Viflux_pm_k04_summer_mean)
print 'variance Kwok04_Jun-Sep_mean:', np.var(Viflux_k04_summer_mean)

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
line3 = slope3 * Vif_s09[mask_s09] + intercept3
print 'PIOMASvsSpreen09:', r_value3, slope3, p_value3
print 'variance PIOMAS:', np.var(Viflux_pm_s09)
print 'variance Spreen09:', np.var(Vif_s09[mask_s09])

ys = np.arange(2003, 2009)
mean_season = np.zeros(len(ys))
mean_season_p = np.zeros(len(ys))
for y in range(len(ys)):
    list_flux_p = []
    list_flux = []
    for i in range(len(dates_pm_s09)):
        date = dates_pm_s09[i]
        if ys[y]==date.year:
            list_flux_p.append(Viflux_pm_s09[i])
            list_flux.append(Vif_s09[i])
    mean_season[y] = np.nanmean(np.array(list_flux))
    mean_season_p[y] = np.nanmean(np.array(list_flux_p))

s, i, r, p, std =  stats.linregress(ys, mean_season)
s_p, i_p, r_p, p_p, std_p =  stats.linregress(ys, mean_season_p)

print s, r**2
print s_p, r_p**2

# PLOTTING
fsize = cm2inch((10,10))
plt.figure(figsize=fsize)

# Kwok04 winter
plt.plot(np.array([-55,650]), np.array([-55,650]), lw = 1.5, c = '#3F5D7D')
plt.plot(Viflux_k04_winter, line1, ls='--', c=tableau20[12], lw=1)
plt.plot(Vif_s09[mask_s09], line3, ls='--', c=tableau20[4], lw=1)
plt.plot(Vif_Cr[mask][:-4], line, ls='--', c=tableau20[2], lw=1)

plt.plot(Viflux_k04_winter, Viflux_pm_k04_winter, 'o', ms=4, markerfacecolor=tableau20[12], markeredgecolor = 'k', markeredgewidth = 0.2, label='Kwok et al. (2004), Oct-May')

# Spreen09
plt.plot(Vif_s09[mask_s09], Viflux_pm_s09, 'o', ms=4, c=tableau20[4], markeredgecolor = 'k', markeredgewidth = 0.2,label='Spreen et al. (2009), Oct-Apr')

# Cryosat
plt.plot(Vif_Cr[mask][:-4], Viflux_pm_cr, 'o', ms=4, c=tableau20[2], markeredgecolor = 'k', markeredgewidth = 0.2,label='Ricker et al. (2017), Oct-Apr')




plt.xlim((-55, 650))
plt.ylim((-55, 650))
plt.title('Monthly', fontsize = 20)
plt.xlabel('Observations', fontsize = '18')
plt.ylabel('PIOMAS' ,fontsize = 18)
plt.xticks( fontsize = 16)
plt.yticks( fontsize = 16)
# plt.title('Monthly ice volume flux, km^3')
# plt.legend(fontsize = 12)
plt.grid( linestyle = '--')
plt.tight_layout()
plt.show()
plt.savefig(OUTDIR+'PIOMAS_Viflux_comparison_scatter.pdf')
