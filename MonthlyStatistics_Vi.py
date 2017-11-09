#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 18:33:30 2017

@author: valeria
"""

from struct import unpack
import numpy as np
import glob, os
import matplotlib.pyplot as plt
import datetime
from scipy import stats


fname_Gr ='./MonthlyIceVolume_PIOMAS_Gr.txt'
fname_Ir = './MonthlyIceVolume_PIOMAS_Ir.txt'

y,m,vi = np.loadtxt(fname_Gr, unpack=True)
vi=vi/1000/1000/1000 #[km^3]
months = range(12)
years = np.arange(1978,2017)


mean_m = np.zeros((12))
max_m = np.zeros((12))
min_m = np.zeros((12))
std_m = np.zeros((12))

for i in months:
    vi_m = []
    for j in range(len(m)):
        if m[j] ==months[i]+1:
            vi_m.append(vi[j])
    vi_m = np.array(vi_m)
    mean_m[i]=vi_m.mean()
    max_m[i]=vi_m.max()
    min_m[i]=vi_m.min()
    std_m[i]=vi_m.std()
    
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, vi_m) 
    print months[i]+1
    print slope1, r_value1, p_value1, std_err1   
    line1 = slope1*years+intercept1
    plt.figure()
    plt.title(str(months[i]+1))
    plt.plot(years,vi_m)
    plt.plot(years,line1)
    plt.show()