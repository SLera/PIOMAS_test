#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 18:04:47 2018

@author: valeria
"""

import numpy as np
from struct import unpack
import ogr, os, sys
import gdal
from gdalconst import *
from osgeo import osr
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime
from scipy import stats
import sys
sys.path.append('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')
import Read_Hi_data

# def create_nc(data):

INDIR_PIOMAS = '/home/valeria/DATA/PIOMAS/v2.1/EASE_grid/rricker'
INDIR_DC_masks = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/DC_regions/DC_regions_masks/'

flist_PIOMAS = [os.path.join(INDIR_PIOMAS, f) for f in os.listdir(INDIR_PIOMAS)]
flist_PIOMAS.sort()

mask_name = 'GreenlandSea_DCmask.npy'
DC_mask = np.load(INDIR_DC_masks+mask_name)
f = open(mask_name[:-4]+'_meanHi.txt', 'w')
for i in range(len(flist_PIOMAS)):
    print i
    date_PIOMAS = Read_Hi_data.extract_date_PIOMAS(flist_PIOMAS[i])
    #прочитать .nc файл с толщинаями PIOMAS
    plat, plon, phi = Read_Hi_data.read_PIOMAS(flist_PIOMAS[i])
    #сделать выборку толщин по маске региона
    phi_reg = Read_Hi_data.extract_region_PIOMAS(phi,DC_mask)
#     #отфильтровать nan
#     phi_reg = phi_reg[~np.isnan(phi_reg)]
#     #найти среднее значение толщины в региона
#     phi_reg_mean = np.mean(phi_reg)
#     #записать строку "дата, средняя толщина в регионе (м)"в текстовый файл
#     line = date_PIOMAS.strftime('%Y%m%d')+'\t'+ str(np.round(phi_reg_mean,2))+ '\n'
#     f.writelines(line)
# f.close()



