#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:11:54 2018

@author: valeria
"""
import numpy as np
from netCDF4 import Dataset

fname = '/home/valeria/DATA/PIOMAS/v2.1/EASE_grid/rricker/piomas_197901.nc'

data = Dataset(fname)