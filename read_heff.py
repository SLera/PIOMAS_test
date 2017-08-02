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
data=np.array(unpacked_bytes).reshape(120,360)


#from scipy.io import FortranFile
#
#ff = FortranFile(FILENAME , 'r')
#a = ff.read_reals(dtype='f8')
#
#
#import numpy as np
#with open(FILENAME,'rb') as nf:
#    for k in xrange(120):
#        data = np.fromfile(nf, dtype=np.float32, count = 120)
