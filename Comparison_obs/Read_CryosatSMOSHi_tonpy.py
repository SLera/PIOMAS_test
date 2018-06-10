#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 13:27:46 2018

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

def gdalwarp (input_file, target_file, epsg, xmin, xmax, ymin, ymax, x_size, y_size):
    print 'gdalwarp -t_srs %s -te %s %s %s %s -tr %s %s -overwrite -of GTiff %s %s' % (epsg, xmin, ymin, xmax, ymax, x_size, y_size, input_file, target_file)
#    ##-999 no data
#    os.system('gdalwarp --config GDAL_DATA "/home/valeria/Programs/miniconda/share/gdal" -r near -t_srs %s -te %s %s %s %s -tr %s %s -overwrite -of GTiff %s %s -wo SAMPLE_GRID=YES -wo SAMPLE_STEPS=1000 -et 0.01 -dstnodata -999' % (epsg, xmin, ymin, xmax, ymax, x_size, y_size, input_file, target_file))
    #nan - no data
    os.system('gdalwarp --config GDAL_DATA "/home/valeria/Programs/miniconda/share/gdal" -r near -t_srs %s -te %s %s %s %s -tr %s %s -overwrite -of GTiff %s %s -et 0.01' % (epsg, xmin, ymin, xmax, ymax, x_size, y_size, input_file, target_file))

def savetonp(fname):
    #tiff to np.array
    dataset = gdal.Open(fname)
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    band = dataset.GetRasterBand(1)
    data = band.ReadAsArray(0, 0, cols, rows)
    data.dump(fname[:-4]+'.npy')
    return data

def create_tif(nparray, output_tiff):
    driver = gdal.GetDriverByName('GTiff')
    outData = driver.Create('temp_5.tif', nparray.shape[1], nparray.shape[0], 1, gdal.GDT_Float32)
    outData.GetRasterBand(1).WriteArray(nparray)
    outData.SetGeoTransform(geotransform_opt)
    outData.SetProjection(wkt_EASE2N)
#    outData.SetProjection(proj4_EASE2N)
    outData.FlushCache()
    del outData
#    gdalwarp('temp_5.tif', output_tiff, target_epsg, target_xmin, target_xmax, target_ymin, target_ymax, target_xsize, target_ysize)
#    savetonp('temp_5.tif')

INDIR = '/home/valeria/DATA/Cryosat2_SMOS/'
fname = INDIR+'cs2smos_ice_thickness_20101115_20101121_v1.3.nc'

##### CONSTS U, V data
xSize = 25000 
ySize = -25000
xCorner = -8987500
yCorner = 8987500
geotransform_opt = [xCorner, xSize, 0, yCorner, 0, ySize]
wkt_EASE2N = 'PROJCS["unnamed",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["latitude_of_center",90],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1],AUTHORITY["epsg","6931"]]'

# # CREATE 25 KM GRID, EASE-2
target_epsg = 'EPSG:6931'
target_xmin = -8987500
target_ymin = -8987500
target_xmax = 8987500
target_ymax = 8987500
target_xsize = 25000
target_ysize = 25000

lat, lon, Hi = read_CryosatSMOS(fname)

create_tif(Hi,'outtif_nan.tif')  

hi = savetonp('temp_5.tif')

plt.figure()
plt.imshow(Hi)
plt.title('Hi')

plt.figure()
plt.imshow(hi)

plt.show()


def tiff_to_latlon(filename, x_im,y_im):
    dataset = gdal.Open(filename, GA_ReadOnly)
    geotransform = dataset.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeigh = geotransform[5]
    
    from osgeo import osr
    # get the existing coordinate system
    old_cs= osr.SpatialReference()
    old_cs.ImportFromWkt(dataset.GetProjectionRef())
    
    # create the new coordinate system
    wgs84_wkt = """
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563,
                AUTHORITY["EPSG","7030"]],
            AUTHORITY["EPSG","6326"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4326"]]"""
    new_cs = osr.SpatialReference()
    new_cs .ImportFromWkt(wgs84_wkt)
    
    # create a transform object to convert between coordinate systems
    transform = osr.CoordinateTransformation(old_cs,new_cs) 
         
    #get the coordinates in lat long
    x = x_im*pixelWidth+originX
    y = y_im*pixelHeigh+originY
    latlong = transform.TransformPoint(x,y) 
    return (latlong[0],latlong[1])

lalo = tiff_to_latlon('temp_5.tif',500,500)

latsEASE = np.zeros((720,720))
lonsEASE = np.zeros((720,720))
for i in range(720):
    print i
    for j in range(720):
        lalo = tiff_to_latlon('temp_5.tif',i,j)
        latsEASE[i,j]= lalo[1]
        lonsEASE[i,j]= lalo[0]

latsEASE.dump('latsEASE2N_25.npy')
lonsEASE.dump('lonsEASE2N_25.npy')