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
INDIR_data = './'



#load regional mask
mask = np.load(INDIR_vars+'Mask_Greenland')

FILENAME = INDIR_data+'heff.H2016'

heffm = read_heff(FILENAME)


