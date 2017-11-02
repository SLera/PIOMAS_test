#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:44:55 2017

@author: valeria
"""

from struct import unpack
import numpy as np
import glob, os
import matplotlib.pyplot as plt

INDIR = './vars/'

#load lon, lat
fname_grid = INDIR+'grid.dat'
a = np.loadtxt(fname_grid)
b = a.flatten()
np.shape(b)
lon = b[0:43200]
lons = lon.reshape(120,360)
lat = b[43200:]
lats = lat.reshape(120,360)

#load kmt
kmt = np.loadtxt(INDIR+'kmt_test.txt')

#region lims [lat0,lat1,lon0,lon1]
REGION = [49,66,300,340]

IrmingerSea_mask = np.zeros((120,360))
IrmingerSea_mask_kmt = np.zeros((120,360))
for i in range(120):
    for j in range(360):
        if (lats[i,j]>REGION[0] and lats[i,j]<REGION[1]):
            if (lons[i,j]>REGION[2] and lons[i,j]<REGION[3]):
                IrmingerSea_mask[i,j] = np.zeros((120,360))[i,j]=1
                if kmt[i,j]>0:
                    IrmingerSea_mask_kmt[i,j]=1
                    
IrmingerSea_mask.dump(INDIR+'LabradorIrminger_mask')
IrmingerSea_mask_kmt.dump(INDIR+'LabradorIrminger_mask')
np.savetxt(INDIR+'LabradorIrminger_mask.txt',IrmingerSea_mask, fmt ='%i')
np.savetxt(INDIR+'LabradorIrminger_mask_kmt.txt',IrmingerSea_mask_kmt, fmt ='%i')

#plot mask    
from mpl_toolkits.basemap import Basemap

plt.figure()
m = Basemap(projection='npstere',boundinglat=48,lon_0=0,resolution='l')
m.drawcoastlines()

(x, y) = m(lons,lats)

m.pcolormesh(x, y, IrmingerSea_mask+IrmingerSea_mask_kmt)
plt.colorbar()
plt.title('Labrador&Irminger')
plt.show()
plt.savefig(INDIR+'Labrador&Irminger Sea mask')