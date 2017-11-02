#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 13:37:05 2017

@author: valeria
"""
from struct import unpack
import numpy as np
import glob, os
import matplotlib.pyplot as plt

#load lon, lat
fname_grid = 'grid.dat'
a = np.loadtxt(fname_grid)
b = a.flatten()
np.shape(b)
lon = b[0:43200]
lons = lon.reshape(120,360)
lat = b[43200:]
lats = lat.reshape(120,360)

#a = heffm[0]
#
#plt.figure()
#plt.imshow(a)
#plt.colorbar()
#plt.show()

#read grid.dat.pop
fname_grid = 'grid.dat.pop'
a = np.loadtxt(fname_grid)
b = a.flatten()
np.shape(b)
#ulon = b[0:(120*360)*1]
#ulons = lon.reshape(120,360)
#ulat = b[(120*360)*1:(120*360)*2]
#ulats = lat.reshape(120,360)
HTN = b[(120*360)*2:(120*360)*3]
HTN = HTN.reshape(120,360)
HTE = b[(120*360)*3:(120*360)*4]
HTE = HTE.reshape(120,360)
#HUS = b[(120*360)*4:(120*360)*5]
#HUS = HUS.reshape(120,360)
#HUW = b[(120*360)*5:(120*360)*6]
#HUW = HUW.reshape(120,360)
#angle= b[(120*360)*6:(120*360)*7]
#angle = angle.reshape(120,360)
c0=0
c1=1.0
p5=0.5
WORK=np.zeros(np.shape(HTN))
WORK[1:,:]=HTN[0:-1,:]
dxt=p5*(HTN+WORK)
dxt[np.where(dxt==c0)]=c1

WORK=np.zeros(np.shape(HTE))
WORK[:,1:]=HTE[:,0:-1]
dyt=p5*(HTE+WORK)
dyt[np.where(dyt==c0)]=c1

dyt.dump('dyt')
dxt.dump('dxt')