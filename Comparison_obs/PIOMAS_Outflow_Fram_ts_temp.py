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
import datetime
from scipy import stats

#import PIOMAS volume flux
INDIR = '/home/valeria/NIERSC/Projects/Bashmachnikov/data/PIOMAS/'
#import ice export from txt
timestamp, Viflux, Sdrift, Hi, Ci = np.loadtxt(INDIR+'PIOMAS_vol_export.txt',skiprows = 1, unpack=True)
days_in_months = [0,31,28,31,30,31,30,31,31,30,31,30,31]

dates = []
years = []
months = []
for i in range(len(timestamp)):
    strng = str(timestamp[i])
    year = int(strng[0:4])
    month = int(strng[4:6])
    Viflux[i]=Viflux[i]*days_in_months[month]
    dates.append(datetime.datetime(year, month, 1))
    years.append(year)
    months.append(month)
    

days_in_months = [0,31,28,31,30,31,30,31,31,30,31,30,31]


years_un = np.arange(1979,2017,1)
months_un = np.arange(1,13,1)


#plt.figure()
#plt.plot(dates,Viflux)
#plt.show()

#import Cryosat volume flux
years_cr = np.arange(2010,2017,1)
#months_cr = np

timestamp_cr, Viflux_cr = np.loadtxt(INDIR+'CRYOSAT_vol_export.txt', unpack=True)

dates_cr =[]
for j in range(len(Viflux_cr)):
    strng = str(timestamp_cr[j])
    year = int(strng[0:4])
    month = int(strng[4:6])
    if Viflux_cr[j]==999:
        Viflux_cr[j]=np.nan
    else:
        Viflux_cr[j]=Viflux_cr[j]*(-1)
    dates_cr.append(datetime.datetime(year, month, 1))
    
#import Kwok_2004 volume flux
timestamp_k04, Viflux_k04 = np.loadtxt(INDIR+'Kwok_vol_export_2004.txt', unpack=True)

dates_k04 =[]
for j in range(len(Viflux_k04)):
    strng = str(timestamp_k04[j])
    year = int(strng[0:4])
    month = int(strng[4:6])
    dates_k04.append(datetime.datetime(year, month, 1))
    
#import Spreen_2009 volume flux
timestamp_s09, Viflux_s09 = np.loadtxt(INDIR+'Spreen_vol_export_2009.txt', unpack=True)

dates_s09 =[]
for j in range(len(Viflux_s09)):
    strng = str(timestamp_s09[j])
    year = int(strng[0:4])
    month = int(strng[4:6])
    if Viflux_s09[j]==-999:
        Viflux_s09[j]=np.nan
    dates_s09.append(datetime.datetime(year, month, 1))

#plot timeseries Cryosat, PIOMAS, Kwock04, Spreen09
ind = np.where(np.array(years)==1991)
i = ind[0][0]

plt.figure(figsize=(15,6))

plt.plot(dates[i:], Viflux[i:], c = 'royalblue', label = 'PIOMAS')
plt.plot(dates_cr[:-4], Viflux_cr[:-4],  c = 'darkorange', label = 'Cryosat-2/SMOS')
plt.plot(dates_k04, Viflux_k04, c = 'plum',  label = 'Kwok2004')
plt.plot(dates_s09, Viflux_s09, c = 'mediumspringgreen',  label = 'Spreen2009')

#plt.plot(dates[i:], Viflux[i:], 'o', ms = 4, markeredgecolor= 'royalblue', markerfacecolor = 'none')
#plt.plot(dates_cr, Viflux_cr,'o', ms = 4,  markeredgecolor = 'darkorange', markerfacecolor = 'none')
#plt.plot(dates_k04, Viflux_k04,'o', ms = 4, markeredgecolor='plum', markerfacecolor = 'none')

plt.ylabel('Monthly ice volume flux, km^3')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
plt.savefig('PIOMAS_Viflux_comparison_ts.pdf')

#create data for scatterplot POIMAS/CRYOSAT
dates_pm_cr = []
Viflux_pm_cr = []
winter_months = [10,11,12,1,2,3,4]

for n in range(0,len(dates_cr[:-4])):
    date_cr = dates_cr[n]
    if date_cr.month in winter_months:
        timelag = dates_cr[n]-np.array(dates)
        timelag_ind = np.where(abs(timelag)==abs(timelag).min())[0][0]
        test = timelag[timelag_ind]
        if test.days !=0:
            print test.days, 'date_cr:', date_cr, 'date_pm:', dates[timelag_ind]
        else:
            dates_pm_cr.append(dates[timelag_ind])
            Viflux_pm_cr.append(Viflux[timelag_ind])
    
mask = ~np.isnan(Viflux_cr)
slope, intercept, r_value, p_value, std_err = stats.linregress(Viflux_cr[mask][:-4],Viflux_pm_cr)     
line = slope*Viflux_cr[mask][:-4]+intercept
print 'PIOMASvsCryosat:',  r_value, slope, p_value
print 'variance PIOMAS:', np.var(Viflux_pm_cr)
print 'variance Cryosat:', np.var(Viflux_cr[mask][:-4])

#create data for scatterplot POIMAS/Kwok2004
dates_pm_k04_winter = []
Viflux_pm_k04_winter = []
dates_k04_winter = []
Viflux_k04_winter = []

dates_pm_k04_summer = []
Viflux_pm_k04_summer = []
dates_k04_summer = []
Viflux_k04_summer = []

winter_months = [10,11,12,1,2,3,4,5]

for n in range(0,len(dates_k04)):
    date_k04= dates_k04[n]
    if date_k04.month in winter_months:
        dates_k04_winter.append(date_k04)
        Viflux_k04_winter.append(Viflux_k04[n])  
        timelag = dates_k04[n]-np.array(dates)
        timelag_ind = np.where(abs(timelag)==abs(timelag).min())[0][0]
        test = timelag[timelag_ind]
        if test.days !=0:
            print test.days, 'date_k04:', date_cr, 'date_pm:', dates[timelag_ind]
        else:
            dates_pm_k04_winter.append(dates[timelag_ind])
            Viflux_pm_k04_winter.append(Viflux[timelag_ind])
    else:
        dates_k04_summer.append(date_k04)
        Viflux_k04_summer.append(Viflux_k04[n])  
        timelag = dates_k04[n]-np.array(dates)
        timelag_ind = np.where(abs(timelag)==abs(timelag).min())[0][0]
        test = timelag[timelag_ind]
        if test.days !=0:
            print test.days, 'date_k04:', date_cr, 'date_pm:', dates[timelag_ind]
        else:
            dates_pm_k04_summer.append(dates[timelag_ind])
            Viflux_pm_k04_summer.append(Viflux[timelag_ind])

slope1, intercept1, r_value1, p_value1, std_er1 = stats.linregress(Viflux_k04_winter,Viflux_pm_k04_winter)   
line1 = slope1*np.array(Viflux_k04_winter)+intercept1
print 'PIOMASvsKwok04:',  r_value1, slope1, p_value1
print 'variance PIOMAS:', np.var(Viflux_pm_k04_winter)
print 'variance Kwok04:', np.var(Viflux_k04_winter)

##find summer_mean
years_k04_summer = np.arange(1992,2000)
Viflux_k04_summer_mean = np.zeros(8)
Viflux_pm_k04_summer_mean = np.zeros(8)

for i in range(len(Viflux_k04_summer)):
    date_k04 = dates_k04_summer[i]
    for j in range(len(years_k04_summer)):
        if date_k04.year == years_k04_summer[j]:
            Viflux_k04_summer_mean[j] += Viflux_k04_summer[i]
            Viflux_pm_k04_summer_mean[j] += Viflux_pm_k04_summer[i]
        
    
slope2, intercept2, r_value2, p_value2, std_er2 = stats.linregress(Viflux_k04_summer_mean,Viflux_pm_k04_summer_mean)   
line2 = slope2*np.array(Viflux_k04_summer_mean)+intercept2
print 'PIOMASvsKwok04_Jun-Sep_mean:',  r_value2, slope2, p_value2
print 'variance PIOMAS_Jun-Sep_mean:', np.var(Viflux_pm_k04_summer_mean)
print 'variance Kwok04_Jun-Sep_mean:', np.var(Viflux_k04_summer_mean)



#create data for scatterplot POIMAS/Spreen09
dates_pm_s09 = []
Viflux_pm_s09 = []
winter_months = [10,11,12,1,2,3,4]

for n in range(0,len(dates_s09)):
    date_s09 = dates_s09[n]
    if date_s09.month in winter_months:
        timelag = dates_s09[n]-np.array(dates)
        timelag_ind = np.where(abs(timelag)==abs(timelag).min())[0][0]
        test = timelag[timelag_ind]
        if test.days !=0:
            print test.days, 'date_s09:', date_s09, 'date_pm:', dates[timelag_ind]
        else:
            dates_pm_s09.append(dates[timelag_ind])
            Viflux_pm_s09.append(Viflux[timelag_ind])
    
mask_s09 = ~np.isnan(Viflux_s09)
slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(Viflux_s09[mask_s09],Viflux_pm_s09)     
line3 = slope3*Viflux_s09[mask_s09]+intercept3
print 'PIOMASvsSpreen09:',  r_value3, slope3, p_value3
print 'variance PIOMAS:', np.var(Viflux_pm_s09)
print 'variance Spreen09:', np.var(Viflux_s09[mask_s09])


#PLOTTING 
plt.figure()
#Cryosat
plt.plot(Viflux_cr[mask][:-4],Viflux_pm_cr, 'o', ms = 5, c = 'darkorange', label = 'Cryosat-2/SMOS, Oct-Apr')
plt.plot(Viflux_cr[mask][:-4],line,ls ='-', c = 'darkorange',lw=0.9)
#Kwok04 winter
plt.plot(Viflux_k04_winter,Viflux_pm_k04_winter, 'o', ms = 5, c = 'plum', label = 'Kwok04, Oct-May')
plt.plot(Viflux_k04_winter,line1,ls ='-', c = 'plum',lw=0.9)
#Kwok04 summer
plt.plot(Viflux_k04_summer_mean,Viflux_pm_k04_summer_mean, '+',ms = 10, c = 'plum', label = 'Kwok04, Jun-Sep, mean')
plt.plot(Viflux_k04_summer_mean,line2,ls ='-', c = 'plum',lw=0.5)
#Spreen09
plt.plot(Viflux_s09[mask_s09],Viflux_pm_s09, 'o', ms = 5, c = 'mediumspringgreen', label = 'Spreen09, Oct-Apr')
plt.plot(Viflux_s09[mask_s09],line3,ls ='-', c = 'mediumspringgreen',lw=0.9)

plt.xlim((-100,600))
plt.ylim((-100,600))
plt.xlabel('Observations')
plt.ylabel('PIOMAS')
plt.title('Monthly ice volume flux, km^3')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
plt.savefig('PIOMAS_Viflux_comparison_scatter.pdf')
