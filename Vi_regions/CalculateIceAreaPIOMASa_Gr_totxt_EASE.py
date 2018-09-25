#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 20:18:52 2017

@author: valeria
"""

from struct import unpack
import numpy as np
import glob, os
import matplotlib.pyplot as plt
import netCDF4
import calendar
import sys
sys.path.append('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')
import Read_Hi_data

def days_in_month(year,month):
    dt = calendar.monthrange(year,month)[1]
    return dt

def calculateVi_PIOMAS(fname):
    year = int(fname[-9:-5])
    month = int(fname[-5:-3])
    data_set = netCDF4.Dataset(fname)
    # date = datetime.date(2010,1,1)+datetime.timedelta(int(time))
    lat = data_set.variables['latitude'][:]
    lon = data_set.variables['longitude'][:]
    Hi = data_set.variables['sit'][:]
    dx = data_set.variables['dx_grid'][:]
    dy = data_set.variables['dy_grid'][:]
    dt = days_in_month(year,month)
    Vi = Hi*dt*dx*dy*1000*1000
    Vi_km = Vi/1000/1000/1000
    data_set.close()
    return lat[144:576,144:576], lon[144:576,144:576], Vi_km[144:576,144:576]

# def extract_region_PIOMAS(data, mask):
#     ind = np.where(mask == 0)
#     data[ind] = np.nan
#     return data

INDIR_data = '/home/valeria/DATA/PIOMAS/v2.1/EASE_grid/rricker/'
OUTDIR = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/output/Greenland/EASE2N_grid/'
# mask = np.load('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test/POIMASCryosat_Greenland/GreenlandSea_mask_Vi_EASE2N.npy')[144:576, 144:576]

months = np.array(['Jan', 'Feb','Mar', 'Apr','May','Jun','Jul','Aug', 'Sep', 'Oct', 'Nov','Dec' ])
months_n = np.arange(1,13)

flist=[]
for root, dirs, ffiles in os.walk(INDIR_data):
  for ffile in ffiles:
      flist.append(os.path.join(root, ffile))
flist.sort()

datestamp = []
Vi_list = []
for f in range(len(flist)):
    fname = flist[f]
    datestamp.append(int(fname[-9:-5]+fname[-5:-3]))
    print datestamp[-1]
    lat, lon, Vip = calculateVi_PIOMAS(fname)
    # Vip_Gr = Read_Hi_data.extract_Greenland(Vip)
    if np.nanmin(Vip)<0:
        print flist[f]
        plt.figure()
        plt.imshow(Vip, vmax = 0.0)
        plt.title(str(datestamp[f]))
        plt.colorbar()
        plt.savefig('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/output/Greenland/Negative_Hi/'+str(datestamp[f])+'_arctic.png')
        plt.close()
#     Vi_list.append(np.nansum(Vip_Gr))
# table = np.column_stack((np.array(datestamp),np.array(Vi_list)))
# np.savetxt(OUTDIR+'PIOMAS_Vi_Greenland_EASE2N.txt', table,  header = 'timestamp, Sea ice volume in km^3', fmt='%1.2f')
#
# plt.figure()
# plt.plot(np.array(Vi_list))
#
# plt.figure()
# plt.imshow(Vip_Gr)
# plt.colorbar()