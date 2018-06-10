#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 13:29:30 2018

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

def read_Cryosat( fname ):
    data_set = netCDF4.Dataset(fname)
    #date = datetime.date(2010,1,1)+datetime.timedelta(int(time))
    lat = data_set.variables['lat'][:]
    lon = data_set.variables['lon'][:]
    Hi = data_set.variables['sea_ice_thickness'][:]
    Hi_unc = data_set.variables['sea_ice_thickness'][:]
    #Hi = np.flipud(Hi)
    data_set.close()
    return  lat, lon, Hi[0]

def read_PIOMAS( fname ):
    data_set = netCDF4.Dataset(fname)
    #date = datetime.date(2010,1,1)+datetime.timedelta(int(time))
    lat = data_set.variables['latitude'][:]
    lon = data_set.variables['longitude'][:]
    Hi = data_set.variables['sit'][:]
    #Hi = np.flipud(Hi)
    data_set.close()
    return  lat, lon, Hi

#INDIR = '/home/valeria/DATA/Cryosat2_SMOS/'
C_fname ='l3c-awi-seaice-cryosat2-ntc-nh25kmEASE2-201201-fv2.0.nc'
P_fname ='piomas_201201.nc'

INDIR = '/home/valeria/DATA/Cryosat2_SMOS/'
fname = INDIR+'cs2smos_ice_thickness_20101115_20101121_v1.3.nc'

plat,plon,Phi = read_PIOMAS(P_fname)
clat,clon,Chi = read_Cryosat(C_fname)

Phi_s = Phi[144:576,144:576]
slat,slon,Shi = read_CryosatSMOS(fname)
d_hi = Chi-Phi_s

GreenlandSea_mask=np.load('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test/POIMASCryosat_Greenland/GreenlandSea_mask_EASE2N.npy')

from mpl_toolkits.basemap import Basemap

#lats_EASE = np.load('latsEASE2N_25.npy')
#lons_EASE = np.load('lonsEASE2N_25.npy')


lats_EASE = clat
lons_EASE = clon
Hi = d_hi
ind = np.where(GreenlandSea_mask[144:576,144:576]==0)
Hi[ind]=np.nan

#m = Basemap(resolution="i",
#            projection='laea', lat_ts=90, lat_0=90., lon_0=0.,
#            llcrnrlon= lons_EASE[576,144], llcrnrlat= lats_EASE[576,144] ,
#            urcrnrlon= lons_EASE[144,576], urcrnrlat= lats_EASE[144,576])

m = Basemap(resolution="i",
            projection='laea', lat_ts=90, lat_0=90., lon_0=0.,
            llcrnrlon= clon[-1,0], llcrnrlat= clat[-1,0] ,
            urcrnrlon= clon[0,-1], urcrnrlat= clat[0,-1])

m.imshow(Hi, origin = 'upper')
m.drawcoastlines()
plt.colorbar()
m.drawmeridians(np.arange(-180,180,10))
m.drawparallels(np.arange(55,80,10))
plt.show()

#m = Basemap(resolution="i",
#            projection='laea', lat_ts=90, lat_0=90., lon_0=0.,
#            llcrnrlon= lons_EASE[200,500], llcrnrlat= lats_EASE[200,500] ,
#            urcrnrlon= lons_EASE[500,200], urcrnrlat= lats_EASE[500,200])
#m.imshow(Hi[200:500,200:500], origin = 'upper')
#m.drawcoastlines()
#plt.colorbar()
#m.drawmeridians(np.arange(-180,180,10))
#m.drawparallels(np.arange(55,80,10))
#plt.show()

#plt.figure()
#m.imshow(GreenlandSea_mask[200:500,200:500],  origin = 'upper')
#m.drawcoastlines()
#plt.show()
#
#plt.figure()
#plt.imshow(lons, origin = 'upper')
#plt.show()