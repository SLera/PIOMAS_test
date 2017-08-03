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
