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

def calc_trend(years, mean_season_Vi):
	##linear regression PIOMAS 2010-2017
	slope_PIOMAS, intercept_PIOMAS, r_PIOMAS, p_PIOMAS, std_PIOMAS =  stats.linregress(years, mean_season_Vi)
	print '1979-2017 PIOMAS trends, slope:', slope_PIOMAS, 'r^2:',r_PIOMAS**2, 'p-value:', p_PIOMAS
	#trend, km^3 per decade 
	Vi_trend_decade = slope_PIOMAS*10
	#trend, % relative to long-term mean per decade 
	long_mean = mean_season_Vi.mean()
	Vi_rel_trend_decade = slope_PIOMAS*10*100/long_mean
	print Vi_trend_decade, '(Trend in Vi flux PIOMAS 1979-2017, km^3 per decade)'
	print Vi_rel_trend_decade, '(Trend in Vi realtive to long-term mean in Vi flux PIOMAS 1979-2017, % per decade)'
	return Vi_trend_decade, Vi_rel_trend_decade

def find_mean_Vi_annual(years, months):
	INDIR = '/home/lera/NIERSC/Projects/Bashmachnikov/data/PIOMAS/'
	# import PIOMAS volume flux
	fname_PIOMAS = INDIR + 'piomas_vol_export_v2.txt'
	fname_PIOMAS = INDIR + 'piomas_vol_export_66.txt'
	Vif_PIOMAS, dates_PIOMAS = V.read_PIOMAS_Vi_flux(fname_PIOMAS)
	mean_season_Vi = np.zeros(len(years))
	for y in range(len(years)):
		list_flux_p = []
		for i in range(len(dates_PIOMAS)):
			date = dates_PIOMAS[i]
			if years[y]==date.year and date.month in months:
				list_flux_p.append(Vif_PIOMAS[i])
		mean_season_Vi[y] = np.nanmean(np.array(list_flux_p))
	return mean_season_Vi

def find_mean_Vi_seasonal(years, months):
	INDIR = '/home/lera/NIERSC/Projects/Bashmachnikov/data/PIOMAS/'
	# import PIOMAS volume flux
	fname_PIOMAS = INDIR + 'piomas_vol_export_v2.txt'
	fname_PIOMAS = INDIR + 'piomas_vol_export_66.txt'
	Vif_PIOMAS, dates_PIOMAS = V.read_PIOMAS_Vi_flux(fname_PIOMAS)
	mean_season_Vi = np.zeros(len(years))
	for y in range(len(years)):
		list_flux_p = []
		for i in range(len(dates_PIOMAS)):
			date = dates_PIOMAS[i]
			if years[y]==date.year and date.month in months[0:-1]:
				list_flux_p.append(Vif_PIOMAS[i])
			if years[y]==(date.year-1) and date.month == months[-1]:
				list_flux_p.append(Vif_PIOMAS[i])
		mean_season_Vi[y] = np.nanmean(np.array(list_flux_p))
	return mean_season_Vi

years = np.arange(1979, 2017)
months_annual =[1,2,3,4,5,6,7,8,9,10,11,12]
months_winter = [1,2,3,4,10,11,12]
#ANNUAL TREND
print 'ANNUAL TREND'
mean_annual_Vi = find_mean_Vi_annual(years,months_annual)
Vi_trend_decade, Vi_rel_trend_decade = calc_trend(years,mean_annual_Vi)
#WINTER TREND
print 'WINTER TREND'
mean_winter_Vi = find_mean_Vi_seasonal(years,months_winter)
Vi_trend_decade, Vi_rel_trend_decade = calc_trend(years,mean_winter_Vi)
