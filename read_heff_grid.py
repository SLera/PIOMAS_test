from struct import unpack
import numpy as np
import glob, os
import matplotlib.pyplot as plt

INDIR = './'

#load lon, lat
fname_grid = INDIR+'grid.dat'
a = np.loadtxt(fname_grid)
b = a.flatten()
np.shape(b)
lon = b[0:43200]
lons = lon.reshape(120,360)
lat = b[43200:]
lats = lat.reshape(120,360)

#load grid mask
kmt = np.loadtxt('kmt_test.txt')

FILENAME = 'heff.H2016'

#bformat= ">%sh"

bformat= "=43200f"

ns = 360*120*4
f = open(FILENAME, 'rb')
byte_arr=f.read(ns)
f.close()
unpacked_bytes = unpack(bformat, byte_arr)
data=np.array(unpacked_bytes).reshape((120,360), order = 'C')

plt.figure()
plt.imshow(data)
plt.show()

#Basemap
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


m = Basemap(projection='npstere',boundinglat=48,lon_0=0,resolution='l')
m.drawcoastlines()

(x, y) = m(lons,lats)

m.pcolormesh(x, y, data)
plt.colorbar()
plt.show()


