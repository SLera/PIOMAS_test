#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 13:27:46 2018

@author: valeria
"""
import netCDF4
import numpy as np
import datetime

def read_PIOMAS_ns( fname ):
    data_set = netCDF4.Dataset(fname)
    time = data_set.variables['time'][:][0]
    date = datetime.date(2010,1,1)+datetime.timedelta(int(time))
    lat = data_set.variables['latitude'][:]
    lon = data_set.variables['longitude'][:]

    u = data_set.variables['ugos'][:].data[0,:,:]
    u_mask = data_set.variables['ugos'][:].mask[0,:,:]
    ind = np.where(u_mask==True)
    u[ind]=np.nan

    v = data_set.variables['vgos'][:].data[0,:,:]
    v_mask = data_set.variables['vgos'][:].mask[0,:,:]
    ind = np.where(v_mask==True)
    v[ind]=np.nan

    return lat, lon, u, v

INDIR = '/home/valeria/DATA/Cryosat2_SMOS/'
fname = INDIR+'cs2smos_ice_thickness_20101115_20101121_v1.3.nc'

data_set = netCDF4.Dataset(fname)

Hi_analysis = data_set.variables['analysis_thickness']
Hi_analysis = data_set.variables['analysis_thickness']