#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 18:04:47 2018

@author: valeria
"""
import PyQt5
import numpy as np
import os
import sys
sys.path.append('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')
import Read_Hi_data

# INDIR_Cr = '/home/valeria/DATA/Cryosat'
# INDIR_PIOMAS = '/home/valeria/DATA/PIOMAS/v2.1/EASE_grid/rricker'

INDIR_Cr = '/home/lera/NIERSC/DATA/Cryosat'
INDIR_PIOMAS = '/home/lera/NIERSC/DATA/PIOMAS/EASE_grid'
flist_Cr = [os.path.join(INDIR_Cr, f) for f in os.listdir(INDIR_Cr)]
flist_Cr.sort()
flist_PIOMAS = [os.path.join(INDIR_PIOMAS, f) for f in os.listdir(INDIR_PIOMAS)]
flist_PIOMAS.sort()
# OUTDIR = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/figs/Greenland/'

months = np.array([10,11,12,1,2,3,4])

chi_months = [[] for x in xrange(len(months))]
chi_u_months = [[] for x in xrange(len(months))]
phi_months = [[] for x in xrange(len(months))]

for i in range(len(flist_PIOMAS)):
    print i
    date_PIOMAS = Read_Hi_data.extract_date_PIOMAS(flist_PIOMAS[i])
    if date_PIOMAS.month in months:
        for j in range(len(flist_Cr)):
            date_Cr = Read_Hi_data.extract_date_Cr(flist_Cr[j])
            if date_Cr == date_PIOMAS:
                clat, clon, chi, chi_u = Read_Hi_data.read_Cryosat(flist_Cr[j])
                plat, plon, phi = Read_Hi_data.read_PIOMAS(flist_PIOMAS[i])
                ind_filter_P = np.isnan(chi)
                phi[ind_filter_P] = np.nan

                # chi = Read_Hi_data.extract_Greenland(chi)
                # chi_u = Read_Hi_data.extract_Greenland(chi_u)
                # phi = Read_Hi_data.extract_Greenland(phi)

                # phi = Read_Hi_data.cut_domain(phi)
                # clat = Read_Hi_data.cut_domain(clat)
                # clon = Read_Hi_data.cut_domain(clon)
                # chi_u = Read_Hi_data.cut_domain(chi_u)
                # chi = Read_Hi_data.cut_domain(chi)

                ind_month = np.where(months == date_PIOMAS.month)[0][0]
                chi_months[ind_month].append(chi)
                chi_u_months[ind_month].append(chi_u)
                phi_months[ind_month].append(phi)

c = np.array(chi_months)
c.dump('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/CryosatHi_full_sorted_monthly')
p = np.array(phi_months)
p.dump('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/PIOMAS_full_sorted_monthly')
cu = np.array(chi_u_months)
cu.dump('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/CryosatUnc_full_sorted_monthly')
clat.dump('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/clat_full')
clon.dump('/home/lera/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/clon_full')