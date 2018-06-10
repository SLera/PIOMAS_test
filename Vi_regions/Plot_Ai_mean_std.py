#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 16:14:12 2017

@author: valeria
"""

from struct import unpack
import numpy as np
import matplotlib.pyplot as plt


from struct import unpack
import numpy as np
import glob, os
import matplotlib.pyplot as plt
import datetime
from scipy import stats


fname_Gr ='./MonthlyIceArea_PIOMAS_Gr.txt'
fname_Ir = './MonthlyIceArea_PIOMAS_Ir.txt'

y,m,vi = np.loadtxt(fname_Gr, unpack=True)
vi=vi/1000/1000/1000 #[km^3]
months = range(12)
years = np.arange(1978,2017)

#mean winter (December-April) Volume
months_winter = np.array([12,1,2,3,4])
v_winter = np.zeros((12, len(years)))
for j in range(len(m)):
    #if (m[j] in months_winter):
    v_winter[int(m[j]-1),int(y[j]-1978)] = vi[j]

plt.figure()    
# row and column sharing
f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter[11]) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1

ax1.plot(years, v_winter[11], 'ro')
ax1.plot(years, v_winter[11], 'r-')
ax1.plot(years, line1, 'r--')
ax1.set_title('Greenland Sea')
s = 'Dec\n'+str(np.round(slope1, decimals = 1))+',p='+str((np.round(p_value1, decimals = 2)))
ax1.text(2000,900, s)

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter[0]) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1


ax2.plot(years, v_winter[0], 'ro')
ax2.plot(years, v_winter[0], 'r-', lw =0.5)
ax2.plot(years, line1, 'r--', lw =0.5)
s = 'Jan\n'+str(np.round(slope1, decimals = 1))+',p='+str((np.round(p_value1, decimals = 2)))
ax2.text(2000,900, s)

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter[1]) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1

ax3.plot(years, v_winter[1], 'ro')
ax3.plot(years, v_winter[1], 'r-', lw =0.5)
ax3.plot(years, line1, 'r--', lw =0.5)
s = 'Feb\n'+str(np.round(slope1, decimals = 1))+',p='+str((np.round(p_value1, decimals = 2)))
ax3.text(2000,900, s)

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter[2]) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1

ax4.plot(years, v_winter[2], 'ro')
ax4.plot(years, v_winter[2], 'r-', lw =0.5)
ax4.plot(years, line1, 'r--', lw =0.5)
s = 'Mar\n'+str(np.round(slope1, decimals = 1))+',p='+str((np.round(p_value1, decimals = 2)))
ax4.text(2000,900, s)

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter[3]) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1

ax5.plot(years, v_winter[3], 'ro')
ax5.plot(years, v_winter[3], 'r-', lw =0.5)
ax5.plot(years, line1, 'r--', lw =0.5)
s = 'Apr'+str(np.round(slope1, decimals = 1))+',p='+str((np.round(p_value1, decimals = 2)))
ax5.text(2000,900, s)

plt.show()

y,m,vi = np.loadtxt(fname_Ir, unpack=True)
vi=vi/1000/1000/1000 #[km^3]
months = range(12)
years = np.arange(1978,2017)

#mean winter (December-April) Volume
months_winter = np.array([12,1,2,3,4])
v_winter = np.zeros((12, len(years)))
for j in range(len(m)):
    #if (m[j] in months_winter):
    v_winter[int(m[j]-1),int(y[j]-1978)] = vi[j]

plt.figure()    
# row and column sharing
f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter[11]) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1

ax1.plot(years, v_winter[11], 'bo')
ax1.plot(years, v_winter[11], 'b-')
ax1.plot(years, line1, 'b--')
ax1.set_title('Irminger/Labrador Seas')
s = 'Dec\n'+str(np.round(slope1, decimals = 1))+',p='+str((np.round(p_value1, decimals = 2)))
ax1.text(2000,400, s)

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter[0]) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1


ax2.plot(years, v_winter[0], 'bo')
ax2.plot(years, v_winter[0], 'b-', lw =0.5)
ax2.plot(years, line1, 'b--', lw =0.5)
s = 'Jan\n'+str(np.round(slope1, decimals = 1))+',p='+str((np.round(p_value1, decimals = 2)))
ax2.text(2000,400, s)

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter[1]) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1

ax3.plot(years, v_winter[1], 'bo')
ax3.plot(years, v_winter[1], 'b-', lw =0.5)
ax3.plot(years, line1, 'b--', lw =0.5)
s = 'Feb\n'+str(np.round(slope1, decimals = 1))+',p='+str((np.round(p_value1, decimals = 2)))
ax3.text(2000,400, s)

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter[2]) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1

ax4.plot(years, v_winter[2], 'bo')
ax4.plot(years, v_winter[2], 'b-', lw =0.5)
ax4.plot(years, line1, 'b--', lw =0.5)
s = 'Mar\n'+str(np.round(slope1, decimals = 1))+',p='+str((np.round(p_value1, decimals = 2)))
ax4.text(2000,400, s)

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years, v_winter[3]) 
print slope1, r_value1, p_value1, std_err1   
line1 = slope1*years+intercept1

ax5.plot(years, v_winter[3], 'bo')
ax5.plot(years, v_winter[3], 'b-', lw =0.5)
ax5.plot(years, line1, 'b--', lw =0.5)
s = 'Apr\n'+str(np.round(slope1, decimals = 1))+',p='+str((np.round(p_value1, decimals = 2)))
ax5.text(2000,400, s)
plt.show()