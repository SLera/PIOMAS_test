#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 18:04:47 2018

@author: valeria
"""

import numpy as np
from netCDF4 import Dataset
import os
import sys
sys.path.append('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')
import Read_Hi_data

def write_nc(OUTDIR,basename, phi_reg, plat, plon):
    ncfile_new = Dataset(OUTDIR+basename+'_PIOMAS_Hi.nc', 'w', format='NETCDF4_CLASSIC')
    dx= ncfile_new.createDimension('dx', len(plon))
    dy = ncfile_new.createDimension('dy', len(plat))
    lon_new = ncfile_new.createVariable('lon', np.float32,('dx','dy'))
    lat_new = ncfile_new.createVariable('lat', np.float32,('dx','dy'))
    str_new = ncfile_new.createVariable('sit', np.float32,('dx', 'dy'))
    #Variable Attributes
    lon_new.units = 'degree_east'
    lon_new.long_name = 'longitude'
    lat_new.units = 'degree_north'
    lat_new.long_name = 'latitude'
    lat_new[:,:] = plat
    lon_new[:,:] = plon
    str_new[:,:] = phi_reg
    ncfile_new.close()
    return

INDIR_PIOMAS = '/home/valeria/DATA/PIOMAS/v2.1/EASE_grid/rricker'
INDIR_DC_masks = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/DC_regions/DC_regions_masks/'
OUTDIR = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/DC_regions/DC_regions_Hi/'

flist_PIOMAS = [os.path.join(INDIR_PIOMAS, f) for f in os.listdir(INDIR_PIOMAS)]
flist_PIOMAS.sort()

mask_name = 'GreenlandSea_DCmask'
DC_mask = np.load(INDIR_DC_masks+mask_name+'.npy')

for i in range(len(flist_PIOMAS)):
    print i
    date_PIOMAS = Read_Hi_data.extract_date_PIOMAS(flist_PIOMAS[i])
    #прочитать .nc файл с толщинаями PIOMAS
    plat, plon, phi = Read_Hi_data.read_PIOMAS(flist_PIOMAS[i])
    #сделать выборку толщин по маске региона
    phi_reg = Read_Hi_data.extract_region_PIOMAS(phi,DC_mask)
    # записать в netscdf
    write_nc(OUTDIR, mask_name, phi_reg, plat, plon)



