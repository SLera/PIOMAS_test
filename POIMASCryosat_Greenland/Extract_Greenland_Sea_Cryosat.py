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


def read_CryosatSMOS( fname ):
    data_set = netCDF4.Dataset(fname)
    #date = datetime.date(2010,1,1)+datetime.timedelta(int(time))
    lat = data_set.variables['latitude'][:]
    lon = data_set.variables['longitude'][:]
    Hi = data_set.variables['analysis_thickness'][:]
    #Hi = np.flipud(Hi)
    data_set.close()
    return  lat, lon, Hi


lats, lons, Hi = read_CryosatSMOS(fname)
REGION = [74,79,340,0]

GreenlandSea_mask = np.zeros((120,360))
GreenlandSea_mask_kmt = np.zeros((120,360))
for i in range(120):
    for j in range(360):
        if (lats[i,j]>REGION[0] and lats[i,j]<REGION[1]):
            if (lons[i,j]>REGION[2] or lons[i,j]<REGION[3]):
                GreenlandSea_mask[i,j]=1
                    
GreenlandSea_mask.dump('GreenlandSea_mask_EASE2N.npy')

