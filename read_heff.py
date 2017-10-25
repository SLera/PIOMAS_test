from struct import unpack
import numpy as np
import glob, os
import matplotlib.pyplot as plt

FILENAME = 'heff.H2016'

#bformat= ">%sh"

bformat= "=43200f"

ns = 360*120*4
f = open(FILENAME, 'rb')
byte_arr=f.read(ns)
f.close()
unpacked_bytes = unpack(bformat, byte_arr)
data=np.array(unpacked_bytes).reshape((360,120), order = 'F')

plt.figure()
plt.imshow(data)
plt.show()


FILENAME = 'iceprod.H2013'

#bformat= ">%sh"

bformat= "=43200f"

ns = 360*120*4
f = open(FILENAME, 'rb')
byte_arr=f.read(ns)
f.close()
unpacked_bytes = unpack(bformat, byte_arr)
data=np.array(unpacked_bytes).reshape((360,120), order = 'F')

data_my=data*60*60*24*31
ind = np.where(data_my==0)
data_my_NaN = data_my
data_my_NaN[ind]=np.nan
plt.figure()
plt.imshow(data_my_NaN, vmin = 2, vmax = 4)
plt.colorbar()
plt.title('Ice productinon, m, Jan 2013')
plt.show()


FILENAME = 'area.H2013'

#bformat= ">%sh"

bformat= "=43200f"

ns = 360*120*4
f = open(FILENAME, 'rb')
byte_arr=f.read(ns)
f.close()
unpacked_bytes = unpack(bformat, byte_arr)
area=np.array(unpacked_bytes).reshape((360,120), order = 'F')


ind = np.where(area<=0)
area_NaN = area
area_NaN[ind]=np.nan
plt.figure()
plt.imshow(area_NaN)
plt.colorbar()
plt.title('Ice area, Jan 2013')
plt.show()
