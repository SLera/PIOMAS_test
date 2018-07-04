#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 13:06:03 2018

@author: valeria
"""
import numpy as np
import datetime
from scipy import stats

def read_PIOMAS_Vi_flux(fname):
    timestamp, Viflux, Sdrift, Hi, Ci = np.loadtxt(fname,skiprows = 1, unpack=True)
    days_in_months = [0,31,28,31,30,31,30,31,31,30,31,30,31]
    dates = []
    for i in range(len(timestamp)):
        strng = str(timestamp[i])
        year = int(strng[0:4])
        month = int(strng[4:6])
        Viflux[i]=Viflux[i]*days_in_months[month]
        dates.append(datetime.datetime(year, month, 1))
    return Viflux, dates

def read_PIOMAS_Vi_flux_coef(fname,c):
    timestamp, Viflux, Sdrift, Hi, Ci = np.loadtxt(fname,skiprows = 1, unpack=True)
    days_in_months = [0,31,28,31,30,31,30,31,31,30,31,30,31]
    dates = []
    for i in range(len(timestamp)):
        strng = str(timestamp[i])
        year = int(strng[0:4])
        month = int(strng[4:6])
        Viflux[i]=Viflux[i]*days_in_months[month]*c
        dates.append(datetime.datetime(year, month, 1))
    return Viflux, dates
    
def read_Cryosat_Vi_flux(fname):
    timestamp_cr, Viflux_cr = np.loadtxt(fname, unpack=True)
    dates_cr =[]
    for j in range(len(Viflux_cr)):
        strng = str(timestamp_cr[j])
        year = int(strng[0:4])
        month = int(strng[4:6])
        if Viflux_cr[j]==999:
            Viflux_cr[j]=np.nan
        else:
            Viflux_cr[j]=Viflux_cr[j]*(-1)
        dates_cr.append(datetime.datetime(year, month, 1))
    return Viflux_cr, dates_cr
    
def read_Kwk_Sprn_Vi_flux(fname):
    timestamp_s09, Viflux_s09 = np.loadtxt(fname, unpack=True)
    dates_s09 =[]
    for j in range(len(Viflux_s09)):
        strng = str(timestamp_s09[j])
        year = int(strng[0:4])
        month = int(strng[4:6])
        if Viflux_s09[j]==-999:
            Viflux_s09[j]=np.nan
        dates_s09.append(datetime.datetime(year, month, 1))
    return Viflux_s09, dates_s09

def sort_by_months(months_cycle, A, dates):
    """Returns values of the parameter A as a list of lists sortes by months 
    -----------
    months_cycle - np.array containing months number for the analysis, e.g. np.array([1,2,3]) for January-March
    A - np.array of a values, e.g. sea ice area
    dates - list of datetime.dates elements corresponding to A-values
    """
    A_montly= []
    months = np.array([d.month for d in dates])
#    years = np.array([d.year for d in dates])
    
    for i in range(len(months_cycle)):
        ind_m = np.where(months==months_cycle[i])
        A_montly.append(A[ind_m])
    return A_montly

def sort_by_years(years_range, A, dates):
    years = np.array([d.year for d in dates])
    ind_y = np.where((years>=years_range[0])&(years<=years_range[-1]))
    A_sorted=A[ind_y]
    dates_sorted = np.array(dates)[ind_y]
    return A_sorted, dates_sorted

def mean_to_txt(months_cycle, A_montly, name):
    mean_cycle = np.zeros(np.shape(months_cycle))
    std = np.zeros(np.shape(months_cycle))
    for i in range(len(months_cycle)):
        mean_cycle[i] = A_montly[i].mean()
        std[i] = A_montly[i].std()
    table = np.column_stack((months_cycle,mean_cycle,std))
    np.savetxt(name, table, fmt='%1.2f')  
    return mean_cycle, std

def calculate_monthly_stats_Vif(years_period, months, Vif, dates, outname):
    #filter by yeras
    Vif_period,dates_period = sort_by_years(years_period, Vif, dates)
    #filter by months
    monthly_Vif = sort_by_months(months,Vif_period,dates_period)
    #calculate mean, std, and write to txt
    mean_monthly_Vif, std_monthly_Vif = mean_to_txt(months, monthly_Vif, outname)
    return monthly_Vif, mean_monthly_Vif, std_monthly_Vif