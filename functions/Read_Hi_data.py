#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 18:04:47 2018

@author: valeria
"""
import netCDF4
import numpy as np
import matplotlib
matplotlib.use('qt5agg')
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib
import datetime

def extract_date_PIOMAS(path):
    return datetime.datetime.strptime(path[-9:-3], '%Y%m')


def extract_date_Cr(path):
    return datetime.datetime.strptime(path[-15:-9], '%Y%m')

def read_Cryosat( fname ):
    data_set = netCDF4.Dataset(fname)
    #date = datetime.date(2010,1,1)+datetime.timedelta(int(time))
    lat = data_set.variables['lat'][:]
    lon = data_set.variables['lon'][:]
    Hi = data_set.variables['sea_ice_thickness'][:][0]
    Hi_unc = data_set.variables['uncertainty'][:][0]
    #Hi = np.flipud(Hi)
    data_set.close()

    return  lat, lon, Hi, Hi_unc

def read_PIOMAS( fname ):
    data_set = netCDF4.Dataset(fname)
    #date = datetime.date(2010,1,1)+datetime.timedelta(int(time))
    lat = data_set.variables['latitude'][:]
    lon = data_set.variables['longitude'][:]
    Hi_eff = data_set.variables['sit'][:]
    sic = data_set.variables['conc'][:]
    sic = sic/100 #fraction instead of %
    #Hi = np.flipud(Hi)
    Hi = Hi_eff/sic
    data_set.close()
    return lat[144:576,144:576], lon[144:576,144:576], Hi[144:576,144:576]

# def calculateVi_PIOMAS( fname ):
#     data_set = netCDF4.Dataset(fname)
#     #date = datetime.date(2010,1,1)+datetime.timedelta(int(time))
#     lat = data_set.variables['latitude'][:]
#     lon = data_set.variables['longitude'][:]
#     Hi = data_set.variables['sit'][:]
#     dx = data_set.variables['dx'][:]
#     dy = data_set.variables['dy'][:]
#     Vi =Hi*dx*dy*t*1000*1000
#     data_set.close()

def extract_Greenland(data):
    GreenlandSea_mask = np.load('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test/POIMASCryosat_Greenland/GreenlandSea_mask_EASE2N.npy')
    GreenlandSea_mask = GreenlandSea_mask[144:576, 144:576]
    ind = np.where(GreenlandSea_mask == 0)
    data[ind] = np.nan
    return data

def extract_Greenland_Vi(data):
    GreenlandSea_mask = np.load('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test/POIMASCryosat_Greenland/GreenlandSea_mask_Vi_EASE2N.npy')
    GreenlandSea_mask = GreenlandSea_mask[144:576, 144:576]
    ind = np.where(GreenlandSea_mask == 0)
    data[ind] = np.nan
    return data

def extract_Irminger_Vi(data):
    IrmingerSea_mask = np.load('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test/POIMASCryosat_Greenland/IrmingerLabrador_mask_Vi_EASE2N.npy')
    IrmingerSea_mask = IrmingerSea_mask[144:576, 144:576]
    ind = np.where(IrmingerSea_mask == 0)
    data[ind] = np.nan
    return data

def extract_region_PIOMAS(data, mask):
    ind = np.where(mask == 0)
    data[ind] = np.nan
    return data

def cut_domain (data):
    return data[242:347,123:260]

def downscale_array(nparray, scale):
    new_size = (np.shape(nparray)[0]/scale,np.shape(nparray)[1]/scale)
    res = np.zeros(new_size)
    res.fill(np.nan)
    for i in range(new_size[0]):
        for j in range(new_size[1]):
            sub = nparray[i*scale:i*scale+scale,j*scale:j*scale+scale]
            res[i,j]= np.nanmean(sub)
    return res

def downscale_array_running_mean(nparray, scale):
    new_size = (np.shape(nparray)[0]-scale+1,np.shape(nparray)[1]-scale+1)
    res = np.zeros(new_size)
    res.fill(np.nan)
    for i in range(new_size[0]):
        for j in range(new_size[1]):
            sub = nparray[i:i+scale,j:j+scale]
            res[i,j]= np.nanmean(sub)
    return res

def plot_map_Greenland(param, limits, clon, clat, outfname, title):
    m = Basemap(resolution="c",
                projection='laea', lat_ts=90, lat_0=90., lon_0=0.,
                llcrnrlon=clon[-1, 0], llcrnrlat=clat[-1, 0],
                urcrnrlon=clon[0, -1], urcrnrlat=clat[0, -1])
    current_cmap = matplotlib.cm.seismic
    current_cmap.set_bad(color='#F1F1F1')
    m.imshow(param, vmin=limits[0], vmax=limits[1], origin='upper', cmap = current_cmap)
    m.drawcoastlines()
    plt.colorbar()
    m.drawmeridians(np.arange(-180, 180, 10), labels=[0, 0, 0, 1], fontsize=12)
    m.drawparallels(np.arange(30, 80, 10), labels=[1, 0, 1, 0], fontsize=12)
    m.fillcontinents(color='#3F5D7D', lake_color='#3F5D7D')
    plt.title(title)
    plt.savefig(outfname)
    plt.close()
    return