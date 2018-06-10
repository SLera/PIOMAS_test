#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 13:24:00 2017

@author: valeria
"""

from struct import unpack
import numpy as np
import glob, os
import matplotlib.pyplot as plt

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

INDIR_vars = './vars/'
INDIR_data = '/home/valeria/DATA/PIOMAS/v2.1/heff/'
OUTDIR = './output/Irminger/'

#load regional mask
mask = np.load(INDIR_vars+'LabradorIrminger_mask')
ind_mask = np.where(mask==0)
kmt = np.loadtxt(INDIR_vars+'kmt_test.txt')
ind_kmt = np.where(kmt<1)

months = np.array(['Jan', 'Feb','Mar', 'Apr','May','Jun','Jul','Aug', 'Sep', 'Oct', 'Nov','Dec' ])
months_n = np.arange(1,13)
#load dy and dy for grid cells (in km)
dxt = np.load(INDIR_vars+'dxt')
dyt = np.load(INDIR_vars+'dyt')

flist=[]
for root, dirs, ffiles in os.walk(INDIR_data):
  for ffile in ffiles:
      flist.append(os.path.join(root, ffile))
flist.sort()
  
plt.figure() 
   
for f in range(len(flist)):
    fname = flist[f]
    year =fname[-4:]
    print year
    heffm = read_heff(fname)
    vi_m = np.zeros(np.shape(months))
    for i in range(len(months)):
        m = months[i]
        vi = heffm[i]*dxt*dyt*1000*1000
        vi[ind_mask]=0
        vi[ind_kmt]=0
        vi=vi/1000/1000/1000 #[km^3]
        #print m, vi.sum()
        vi_m[i]=vi.sum()
    #    plt.figure()
    #    plt.imshow(vi,vmin=0.1)
    #    plt.colorbar()
    #    plt.title('Vi,km^3'+m+'max = '+str(vi.max()) )
    #    plt.show()
    
    np.savetxt(OUTDIR+'LabIrmingerSea'+year+'.txt', vi_m, header = 'Monthly (Jan-Dec) Sea ice volume in km^3',fmt='%1.9f')
    
    #plt.figure()
    plt.plot(vi_m)
plt.xticks(np.arange(len(months)),months)
plt.ylabel('Ice volume, km^3')
plt.grid()
plt.title('1978-2016')
plt.show()