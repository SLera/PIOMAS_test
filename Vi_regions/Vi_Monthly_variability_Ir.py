#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 16:41:38 2017

@author: valeria
"""

from struct import unpack
import numpy as np
import glob, os
import matplotlib.pyplot as plt

INDIR = './output/Irminger/'

flist = os.listdir(INDIR)
flist.sort()

months = np.array(['Jan', 'Feb','Mar', 'Apr','May','Jun','Jul','Aug', 'Sep', 'Oct', 'Nov','Dec' ])
months_n = np.arange(1,13)

for m in range(len(months_n)):
    monthly_ts = np.zeros((len(flist)))
    years = np.zeros((len(flist)))
    for f in range(len(flist)):
        year = flist[f][-8:-4]
        vi=np.loadtxt(INDIR+flist[f])
        monthly_ts[f]=vi[m]
        years[f]=int(year)
    mean = monthly_ts.mean()
    rel_std = (monthly_ts.std()*100)/mean
    plt.figure()
    plt.title('1978-2016, '+months[m])
    plt.plot(years,monthly_ts)
    plt.ylim((-10,500))
    plt.grid()
    plt.ylabel('Ice volume, km^3')
    plt.text(1978,400,('%.2f'%round(mean,2))+'+-'+'2x'+('%.1f'%round(rel_std,1))+'%')
    plt.show()
    plt.savefig(str(m+1)+'_'+'1978-2016 Lab&Irm')