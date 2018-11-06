#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 16:14:12 2017

@author: valeria
"""

from struct import unpack
import numpy as np
import matplotlib
matplotlib.use('qt5agg')
import matplotlib.pyplot as plt
from struct import unpack
import numpy as np
import glob, os
import matplotlib.pyplot as plt
import datetime
from scipy import stats


fname_Gr ='/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/MonthlyIceVolume_PIOMAS_Gr.txt'
fname_Ir = '/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/MonthlyIceVolume_PIOMAS_Ir.txt'

fname_Gr = '/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/output/Greenland/EASE2N_grid_v2/PIOMAS_Vi_Greenland_EASE2N.txt'
fname_Ir = '/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/output/Greenland/EASE2N_grid_v2/PIOMAS_Vi_IrmingerLabrador_EASE2N.txt'


y,m,vi = np.loadtxt(fname_Gr, unpack=True)
months = range(12)
years = np.arange(1979,2017)

#mean winter (December-April) Volume
months_winter = np.array([12,1,2,3,4])
v_winter = np.zeros(len(years))
for i in range(len(years)):
    #print years[i]
    winter_year = 0
    for j in range(len(m)):
        if (y[j] == years[i]) and (m[j] in months_winter[1:]):
            winter_year += vi[j]
        if (y[j] == years[i] - 1) and (m[j] == 12):
            winter_year += vi[j]

    if i == 0:
        v_winter[i] = winter_year / 4
    else:
        v_winter[i] = winter_year / 5

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1
    
fig, ax1 = plt.subplots()
ax1.plot(years, v_winter, 'ro', ms=5, label = 'Greenland Sea')
ax1.plot(years, v_winter, 'r-')
ax1.plot(years, line1, 'r--', lw =2)
ax1.set_ylabel('Greenaland Sea winter sea ice volume, km^3', color='r')
#ax1.yaxis.set_ticks(np.arange(500, 1501, 250))
ax1.tick_params('y', colors='r')
#plt.xticks(np.arange(len(months_name)),months_name)
plt.legend(loc=2)
plt.title('Trends in winter (Dec-Apr) sea ice volume')

#Irminger Sea
y,m,vi = np.loadtxt(fname_Ir, unpack=True)
months = range(12)
years = np.arange(1979,2017)

#mean winter (December-April) Volume
months_winter = np.array([12,1,2,3,4])
v_winter = np.zeros(len(years))
for i in range(len(years)):
    #print years[i]
    winter_year = 0
    for j in range(len(m)):
        if (y[j]==years[i]) and (m[j] in months_winter[1:]):
            winter_year += vi[j]
        if (y[j]==years[i]-1) and (m[j] == 12):
            winter_year += vi[j]

    if i == 0:
        v_winter[i] = winter_year/4
    else:
        v_winter[i] = winter_year / 5

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1
    
ax2 = ax1.twinx()
ax2.plot(years, v_winter, 'bo', ms=5, label = 'Irminger/Labrador Seas')
ax2.plot(years, v_winter, 'b-')
ax2.plot(years, line1, 'b--', lw = 1)
ax2.set_ylabel('Irminger/Labrador Seas winter sea ice volume, km^3', color='b')
ax2.tick_params('y', colors='b')
#ax2.yaxis.set_ticks(np.arange(0, 101, 25))
#plt.show()
plt.savefig('ViWinterTrends_EASE.pdf', dpi=400)
