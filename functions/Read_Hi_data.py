#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 18:04:47 2018

@author: valeria
"""
import netCDF4
import numpy as np

def read_Cryosat( fname ):
    data_set = netCDF4.Dataset(fname)
    #date = datetime.date(2010,1,1)+datetime.timedelta(int(time))
    lat = data_set.variables['lat'][:]
    lon = data_set.variables['lon'][:]
    Hi = data_set.variables['sea_ice_thickness'][:][0]
    Hi_unc = data_set.variables['uncertainty'][:]
    #Hi = np.flipud(Hi)
    data_set.close()

    return  lat, lon, Hi, Hi_unc

def read_PIOMAS( fname ):
    data_set = netCDF4.Dataset(fname)
    #date = datetime.date(2010,1,1)+datetime.timedelta(int(time))
    lat = data_set.variables['latitude'][:]
    lon = data_set.variables['longitude'][:]
    Hi = data_set.variables['sit'][:]
    #Hi = np.flipud(Hi)
    data_set.close()
    return lat[144:576,144:576], lon[144:576,144:576], Hi[144:576,144:576]

def extract_Greenland(data):
    GreenlandSea_mask = np.load('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test/POIMASCryosat_Greenland/GreenlandSea_mask_EASE2N.npy')
    GreenlandSea_mask = GreenlandSea_mask[144:576, 144:576]
    ind = np.where(GreenlandSea_mask == 0)
    data[ind] = np.nan
    return data

def cut_domain (data):
    return data[242:347,123:260]
