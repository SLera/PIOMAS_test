#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 15:59:47 2018

@author: valeria
"""
import netCDF4
import numpy as np
import math
import os
import gdal 
from gdalconst import * 
from osgeo import osr
from matplotlib import pyplot as plt

def read_CryosatSMOS( fname ):
    data_set = netCDF4.Dataset(fname)
    #date = datetime.date(2010,1,1)+datetime.timedelta(int(time))
    lat = data_set.variables['latitude'][:]
    lon = data_set.variables['longitude'][:]
    Hi = data_set.variables['analysis_thickness'][:]
    #Hi = np.flipud(Hi)
    data_set.close()
    return  lat, lon, Hi

INDIR = '/home/valeria/DATA/Cryosat2_SMOS/'
fname = INDIR+'cs2smos_ice_thickness_20101115_20101121_v1.3.nc'

lats, lons, Hi = read_CryosatSMOS(fname)
REGION = [60, 82,-45, 20]

GreenlandSea_mask = np.zeros((720,720))

for i in range(720):
    print i
    for j in range(720):
        if (lats[i,j]>REGION[0] and lats[i,j]<REGION[1]):
            if ((lons[i,j]>REGION[2] and lons[i,j]<0) or (lons[i,j]<REGION[3] and lons[i,j]>0)):
                print lons[i,j]
                GreenlandSea_mask[i,j]=1
                print GreenlandSea_mask[i,j]
                    
GreenlandSea_mask.dump('GreenlandSea_mask_EASE2N.npy')

from mpl_toolkits.basemap import Basemap

#lats_EASE = np.load('latsEASE2N_25.npy')
#lons_EASE = np.load('lonsEASE2N_25.npy')


lats_EASE = lats
lons_EASE = lons

m = Basemap(resolution="i",
            projection='laea', lat_ts=90, lat_0=90., lon_0=0.,
            llcrnrlon= lons_EASE[500,200], llcrnrlat= lats_EASE[500,200] ,
            urcrnrlon= lons_EASE[200,500], urcrnrlat= lats_EASE[200,500])
m.imshow(Hi[200:500,200:500], origin = 'upper')
m.drawcoastlines()
plt.show()

plt.figure()
m.imshow(GreenlandSea_mask[200:500,200:500],  origin = 'upper')
m.drawcoastlines()
plt.show()

plt.figure()
plt.imshow(lons, origin = 'upper')
plt.show()