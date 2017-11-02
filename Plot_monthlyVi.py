#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 14:21:14 2017

@author: valeria
"""

from struct import unpack
import numpy as np
import glob, os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def read_heff(FILENAME):
    bformat= "=518400f"
    months = np.arange(1,13)
    ns = 360*120*4*len(months)
    f = open(FILENAME, 'rb')
    byte_arr=f.read(ns)
    f.close()
    unpacked_bytes = unpack(bformat, byte_arr)
    heffm = []
    for i in range(len(months)):
        ind1=(120*360)*i
        ind2=(120*360)*(i+1)
        #data=np.array(unpacked_bytes[ind1:ind2]).reshape((120,360), order = 'C')
        data=np.array(unpacked_bytes[ind1:ind2])
        heffm.append(data.reshape((120,360), order = 'C'))
    return heffm

def plot_POIMAS_map(PIOMAS_data,title):
    #load lon, lat
    fname_grid = './vars/grid.dat'
    a = np.loadtxt(fname_grid)
    b = a.flatten()
    np.shape(b)
    lon = b[0:43200]
    lons = lon.reshape(120,360)
    lat = b[43200:]
    lats = lat.reshape(120,360)
    plt.figure()
    m = Basemap(projection='npstere',boundinglat=48,lon_0=0,resolution='l')
    m.drawcoastlines()
    (x, y) = m(lons,lats)
    m.pcolormesh(x, y, PIOMAS_data)
    plt.colorbar()
    plt.title(title)
    plt.show()
    return 'plot'

INDIR_vars = './vars/'
INDIR_data = './'

months = np.array([1,2,3,4,5,6,7,8,9,10,11,12])

#load regional mask
mask = np.load(INDIR_vars+'Mask_Greenland')
#load dy and dy for grid cells (in km)
dxt = np.load(INDIR_vars+'dxt')
dyt = np.load(INDIR_vars+'dyt')

fname = INDIR_data+'heff.H2016'
year =fname[-4:]
heffm = read_heff(fname)

for i in range(len(months)):
    m = months[i]
    vi = heffm[i]*dxt*dyt*1000*1000
    plot_POIMAS_map(vi,str(months[i]))
