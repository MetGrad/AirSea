#!/usr/bin/env python
# coding: utf-8

# In[91]:

import sys 
import math
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import scipy as sp
from scipy import stats
import cartopy.feature as cfeature
import metpy.calc as mpcalc
import metpy
from metpy.units import units
import pandas as pd 
import statistics
import traceback

#Plots: ocean currents, Wind speed, wind stress curl, Surface vort/curl, ocean divergence


# In[92]:
fp = 'D:/Tutor/FSU_Anna/'

#Variables
#March 26, 2018

#latitude, longitude
#mean_observation_time

#wind_speed in m/s, wind_dir in -180 to 180, 10m equiv neutral
#u_current, v_current in m/s 
#u_current_std, v_current_std in m/s 

#tau_x, tau_y, description m/s, units Pa - downward stress at sea sfc 

#'du_current_dx', 'du_current_dy', 'dv_current_dx', 'dv_current_dy' in 1/s - 
#East component of surface current derivative (dus), North component of surface current derivative (dvs)
#'du_current_dx_std', 'du_current_dy_std', 'dv_current_dx_std', 'dv_current_dy_std' in 1/s 

#'surface_current_relative_vorticity' in 1/s 
#'surface_current_relative_vorticity_std' in 1/s 
#'surface_current_divergence' in 1/s 
#'surface_current_divergence_std' in 1/s 

#'surface_current_strain_rate' in 1/s 
#'coriolis_parameter' in 1/s 

# 'wind_stress_curl' in  N m-3
#'dtau_x_dx', 'dtau_x_dy', 'dtau_y_dx', 'dtau_y_dy' in  N m-3 - 
#East component of wind stress derivative in the e for dTaux, North component of wind stress derivative in the in dTauY

ds4 = xr.open_dataset(f'{fp}20180326_084656_0238-0296_line04.L2.nc')
ds7 = xr.open_dataset(f'{fp}20180326_084656_0377-0419_line07.L2.nc')

Lat4 = ds4.latitude.values
Lon4 = ds4.longitude.values

Lat7 = ds7.latitude.values
Lon7 = ds7.longitude.values

uCur4 = ds4.u_current.values
vCur4 = ds4.v_current.values

uCur7 = ds7.u_current.values
vCur7 = ds7.v_current.values

curMag4 = np.sqrt(uCur4**2 + vCur4**2)
curMag7 = np.sqrt(uCur7**2 + vCur7**2)

numX4 = len(Lon4)
numY4 = len(Lat4)

numX7 = len(Lon7)
numY7 = len(Lat7)

# 'wind_stress_curl' in  N m-3
#'surface_current_divergence' in 1/s 
#'surface_current_relative_vorticity' in 1/s 
#wind_speed in m/s, wind_dir in -180 to 180, 10m equiv neutral

windSp4 = ds4.wind_speed.values 
windDir4 = ds4.wind_dir.values

windSp7 = ds7.wind_speed.values 
windDir7 = ds7.wind_dir.values

wDirMath4 = 90 - windDir4 * (np.pi/180)
wSpeed4 = np.sqrt((windSp4*np.cos(wDirMath4*(180/np.pi)))**2 + (windSp4*np.sin(wDirMath4*(180/np.pi)))**2)

wDirMath7 = 90 - windDir7 * (np.pi/180)
wSpeed7 = np.sqrt((windSp7*np.cos(wDirMath7*(180/np.pi)))**2 + (windSp7*np.sin(wDirMath7*(180/np.pi)))**2)

windU4 = windSp4 * np.cos(wDirMath4 * (180/np.pi))
windV4 = windSp4 * np.sin(wDirMath4 * (180/np.pi))

windU7 = windSp7 * np.cos(wDirMath7 * (180/np.pi))
windV7 = windSp7 * np.sin(wDirMath7 * (180/np.pi))

U4_vec = windU4/windSp4
V4_vec = windV4/windSp4

U7_vec = windU7/windSp7
V7_vec = windV7/windSp7

U4_perp = -V4_vec
V4_perp = U4_vec

U7_perp = -V7_vec
V7_perp = U7_vec

perpCurrent4 = U4_perp*uCur4 + V4_perp*vCur4
perpCurrent7 = U7_perp*uCur7 + V7_perp*vCur7

# 'wind_stress_curl' in  N m-3
#'surface_current_divergence' in 1/s 
#'surface_current_relative_vorticity' in 1/s


stressCurl4 = ds4.wind_stress_curl.values
stressCurl7 = ds7.wind_stress_curl.values

#tau_x, tau_y, description m/s, units Pa - downward stress at sea sfc 


tauX4 = ds4.tau_x.values
tauY4 = ds4.tau_y.values

tauX7 = ds7.tau_x.values
tauY7 = ds7.tau_y.values

stressMag4 = np.sqrt(tauX4**2 + tauY4**2)
stressMag7 = np.sqrt(tauX7**2 + tauY7**2)

  

#3rd dim: dif in y [0] or dif in x [0]

wind67_curdif4 = np.zeros((numY4, numX4, 2)) * np.nan
wind67_curdifperp4 = np.zeros((numY4, numX4, 2)) * np.nan
wind67_stressmagdif4 = np.zeros((numY4, numX4, 2)) * np.nan
wind67_stresscurldif4 = np.zeros((numY4, numX4, 2)) * np.nan

wind67_curdif7 = np.zeros((numY7, numX7, 2)) * np.nan
wind67_curdifperp7 = np.zeros((numY7, numX7, 2)) * np.nan
wind67_stressmagdif7 = np.zeros((numY7, numX7, 2)) * np.nan
wind67_stresscurldif7 = np.zeros((numY7, numX7, 2)) * np.nan

#wind910_curdif = np.zeros((numY4, numX4, 2)) * np.nan
#wind910_curdifperp = np.zeros((numY4, numX4, 2)) * np.nan
#wind910_stressmagdif = np.zeros((numY4, numX4, 2)) * np.nan
#wind910_stresscurldif = np.zeros((numY4, numX4, 2)) * np.nan


#do wind bins for ds4
#WIND SPEED: high = 9 - 10m/s, med = 9 - 7 m/s

# Update for loop to use current magnitude instead of components
#   curMag4=np.sqrt(uCur4**2+vCur4**2)
# Update for loop to use stress magnitude
#   stressMag4=np.sqrt(tauX4**2+tauY4**2)
cur_threshhold = 0.1

for x in range(0, len(Lon4) - 1):
    for y in range(0, len(Lat4) - 1):
        if np.isnan(wSpeed4[y, x]) == False:
            magCurrent_x = curMag4[y, x+1] - curMag4[y, x]
            magCurrent_y = curMag4[y+1, x] - curMag4[y, x]
            
            stressmag_difx = stressMag4[y, x+1] - stressMag4[y, x]
            stressmag_dify = stressMag4[y+1, x] - stressMag4[y, x]
            
            stresscurl_difx = stressCurl4[y, x+1] - stressCurl4[y, x]
            stresscurl_dify = stressCurl4[y+1, x] - stressCurl4[y, x]
            
            perpCurrent_x = perpCurrent4[y, x+1] - perpCurrent4[y, x]
            perpCurrent_y = perpCurrent4[y+1, x] - perpCurrent4[y, x]
                    
            if wSpeed4[y, x] > 6.0 and wSpeed4[y, x] <= 7.0:                   
                wind67_curdif4[y, x, 0] = magCurrent_y                          
                wind67_curdif4[y, x, 1] = magCurrent_x
                wind67_curdifperp4[y, x, 0] = perpCurrent_y
                wind67_curdifperp4[y, x, 1] = perpCurrent_x
                wind67_stressmagdif4[y, x, 0] = stressmag_dify
                wind67_stressmagdif4[y, x, 1] = stressmag_difx
                wind67_stresscurldif4[y, x, 0] = stresscurl_dify
                wind67_stresscurldif4[y, x, 1] = stresscurl_difx
                
            else:
                continue
            
'''
elif wSpeed4[y, x] > 9.0 and wSpeed4[y, x] <= 10.0:
    wind910_curdif[y, x, 0] = magCurrent_y
    wind910_curdif[y, x, 1] = magCurrent_x
    wind910_curdifperp[y, x, 0] = perpCurrent_y
    wind910_curdifperp[y, x, 1] = perpCurrent_x
    wind910_stressmagdif[y, x, 0] = stressmag_dify
    wind910_stressmagdif[y, x, 1] = stressmag_difx
    wind910_stresscurldif[y, x, 0] = stresscurl_dify
    wind910_stresscurldif[y, x, 1] = stresscurl_difx
'''

for x in range(0, len(Lon7) - 1):
    for y in range(0, len(Lat7) - 1):
        if np.isnan(wSpeed7[y, x]) == False:
            magCurrent_x = curMag7[y, x+1] - curMag7[y, x]
            magCurrent_y = curMag7[y+1, x] - curMag7[y, x]
            
            stressmag_difx = stressMag7[y, x+1] - stressMag7[y, x]
            stressmag_dify = stressMag7[y+1, x] - stressMag7[y, x]
            
            stresscurl_difx = stressCurl7[y, x+1] - stressCurl7[y, x]
            stresscurl_dify = stressCurl7[y+1, x] - stressCurl7[y, x]
            
            perpCurrent_x = perpCurrent7[y, x+1] - perpCurrent7[y, x]
            perpCurrent_y = perpCurrent7[y+1, x] - perpCurrent7[y, x]
                    
            if wSpeed7[y, x] > 6.0 and wSpeed7[y, x] <= 7.0:                   
                wind67_curdif7[y, x, 0] = magCurrent_y                          
                wind67_curdif7[y, x, 1] = magCurrent_x
                wind67_curdifperp7[y, x, 0] = perpCurrent_y
                wind67_curdifperp7[y, x, 1] = perpCurrent_x
                wind67_stressmagdif7[y, x, 0] = stressmag_dify
                wind67_stressmagdif7[y, x, 1] = stressmag_difx
                wind67_stresscurldif7[y, x, 0] = stresscurl_dify
                wind67_stresscurldif7[y, x, 1] = stresscurl_difx
                
            else:
                continue



### Flatten
#wind67_curdifx4 = wind67_curdif4[:, :, 1].flatten()
#wind67_curdify4 = wind67_curdif4[:, :, 0].flatten()
wind67_curdifperpx4 = wind67_curdifperp4[:, :, 1].flatten()
wind67_curdifperpy4 = wind67_curdifperp4[:, :, 0].flatten()
#wind67_stressmagdifx4 = wind67_stressmagdif4[:, :, 1].flatten()
#wind67_stressmagdify4 = wind67_stressmagdif4[:, :, 0].flatten()
wind67_stresscurldifx4 = wind67_stresscurldif4[:, :, 1].flatten()
wind67_stresscurldify4 = wind67_stresscurldif4[:, :, 0].flatten()


wind67_curdifperpx7 = wind67_curdifperp7[:, :, 1].flatten()
wind67_curdifperpy7 = wind67_curdifperp7[:, :, 0].flatten()
wind67_stresscurldifx7 = wind67_stresscurldif7[:, :, 1].flatten()
wind67_stresscurldify7 = wind67_stresscurldif7[:, :, 0].flatten()

'''
wind910_curdifx = wind910_curdif[:, :, 1].flatten()
wind910_curdify = wind910_curdif[:, :, 0].flatten()
wind910_curdifperpx = wind910_curdifperp[:, :, 1].flatten()
wind910_curdifperpy = wind910_curdifperp[:, :, 0].flatten()
wind910_stressmagdifx = wind910_stressmagdif[:, :, 1].flatten()
wind910_stressmagdify = wind910_stressmagdif[:, :, 0].flatten()
wind910_stresscurldifx = wind910_stresscurldif[:, :, 1].flatten()
wind910_stresscurldify = wind910_stresscurldif[:, :, 0].flatten()
'''

x67_4 = wind67_curdifperpx4
y67_4 = wind67_stresscurldifx4

x67_7 = wind67_curdifperpx7
y67_7 = wind67_stresscurldifx7

'''
x910 = wind910_curdifperpx
y910 = wind910_stresscurldifx
'''
'''
mask = ~np.isnan(wind78_curdifx) & ~np.isnan(wind78_stresscurldifx)
slope, intercept, rvalue, pvalue, stderr = stats.linregress(wind78_curdifx[mask], wind78_stresscurldifx[mask])
res = stats.linregress(wind78_curdifx[mask], wind78_stresscurldifx[mask])
'''

'''
figure(figsize = (10, 5), dpi = 300)
### Scatter Plot
plt.scatter(x67_7, y67_7, color = 'tab:blue', alpha = 0.5, label = 'Wind 6-7 m/s')
#plt.scatter(x910, y910, color = 'tab:red', alpha = 0.5, label = 'wind 9-10 m/s')
plt.legend(loc = "upper right")
#plt.plot(wind78_curdify[mask], res.intercept + res.slope*wind78_curdify[mask], 'g', alpha = 1)
plt.xlabel('Difference in current (m/s)')
plt.ylabel('Difference in stress Curl (N/m^3)')
plt.show()
#plt.savefig(f'{fp}test_high.png', dpi = 300)
'''


############################### START HERE ####################################
###### Create new arrays with nans removed to reduce loop time
# You will have to do this for each wind bin (6-7, 9-10) for each dataset, for x and y,
# for the stress curl, stress magnitude, and current difference
wind67_curdifperpx4 = wind67_curdifperpx4[~(np.isnan(wind67_curdifperpx4))]
wind67_curdifperpy4 = wind67_curdifperpy4[~(np.isnan(wind67_curdifperpy4))]
wind67_stresscurldifx4 = wind67_stresscurldifx4[~(np.isnan(wind67_stresscurldifx4))]
wind67_stresscurldify4 = wind67_stresscurldify4[~(np.isnan(wind67_stresscurldify4))]

wind67_curdifperpx7 = wind67_curdifperpx7[~(np.isnan(wind67_curdifperpx7))]
wind67_curdifperpy7 = wind67_curdifperpy7[~(np.isnan(wind67_curdifperpy7))]
wind67_stresscurldifx7 = wind67_stresscurldifx7[~(np.isnan(wind67_stresscurldifx7))]
wind67_stresscurldify7 = wind67_stresscurldify7[~(np.isnan(wind67_stresscurldify7))]

bin_n35_n25_x_67 = 0
bin_n25_n15_x_67 = 0
bin_n15_n05_x_67 = 0
bin_n05_00_x_67 = 0
bin_00_05_x_67 = 0
bin_05_15_x_67 = 0
bin_15_25_x_67 = 0
bin_25_35_x_67 = 0


# 287, 242
# Create a meshgrid, and use a simpler variable name to refer to the current difference
# and stress magnitude/stress curl. MAKE SURE TO LIST THE CURRENT DIFFERENCE 
# SECOND, and the output variables (y4, x4) should be named in order y, x
y4, x4 = np.meshgrid(wind67_stresscurldifx4, wind67_curdifperpx4)
#x4 = wind67_curdifperpx4
#y4 = wind67_stresscurldifx4
# testy changes in the 2nd dim (axis 1)
# testx changes in the 1st dim (axis 0)

for x in range(len(wind67_curdifperpx4)):
    
    # Now, use the current difference array you just made
    if x4[x, 0] < -2.5:
        bin_n35_n25_x_67 += 1

    elif x4[x, 0] >= -2.5 and x4[x, 0] < -1.5:            
    #elif wind67_curdifperp[x] >= -2.5 and wind67_curdifperp[x] < -1.5:
        bin_n25_n15_x_67 += 1
        
    elif x4[x, 0] >= -1.5 and x4[x, 0] < -0.5:
    #elif wind67_curdifperp[x] >= -1.5 and wind67_curdifperp[x] < -0.5:
        bin_n15_n05_x_67 += 1
        
    elif x4[x, 0] >= -0.5 and x4[x, 0] < 0.0:
    #elif wind67_curdifperp[x] >= -0.5 and wind67_curdifperp[x] < 0.5:
        bin_n05_00_x_67 += 1
        
    elif x4[x, 0] >= 0.0 and x4[x, 0] < 0.5:
    #elif wind67_curdifperp[x] >= -0.5 and wind67_curdifperp[x] < 0.5:
        bin_00_05_x_67 += 1
        
    elif x4[x, 0] >= 0.5 and x4[x, 0] < 1.5:
    #elif wind67_curdifperp[x] >= 0.5 and wind67_curdifperp[x] < 1.5:
        bin_05_15_x_67 += 1
        
    elif x4[x, 0] >= 1.5 and x4[x, 0] < 2.5:
    #elif wind67_curdifperp[x] >= 1.5 and wind67_curdifperp[x] < 2.5:
        bin_15_25_x_67 += 1
        
    else:
        bin_25_35_x_67 += 1
    

y7, x7 = np.meshgrid(wind67_stresscurldifx7, wind67_curdifperpx7)
# testy changes in the 2nd dim (axis 1)
# testx changes in the 1st dim (axis 0)

for x in range(len(wind67_curdifperpx7)):
    
    if x7[x, 0] < -2.5:
        bin_n35_n25_x_67 += 1

    elif x7[x, 0] >= -2.5 and x7[x, 0] < -1.5:            
    #elif wind67_curdifperp[x] >= -2.5 and wind67_curdifperp[x] < -1.5:
        bin_n25_n15_x_67 += 1
        
    elif x7[x, 0] >= -1.5 and x7[x, 0] < -0.5:
    #elif wind67_curdifperp[x] >= -1.5 and wind67_curdifperp[x] < -0.5:
        bin_n15_n05_x_67 += 1
        
    elif x7[x, 0] >= -0.5 and x7[x, 0] < 0.0:
    #elif wind67_curdifperp[x] >= -0.5 and wind67_curdifperp[x] < 0.5:
        bin_n05_00_x_67 += 1
        
    elif x7[x, 0] >= 0.0 and x7[x, 0] < 0.5:
    #elif wind67_curdifperp[x] >= -0.5 and wind67_curdifperp[x] < 0.5:
        bin_00_05_x_67 += 1
        
    elif x7[x, 0] >= 0.5 and x7[x, 0] < 1.5:
    #elif wind67_curdifperp[x] >= 0.5 and wind67_curdifperp[x] < 1.5:
        bin_05_15_x_67 += 1
        
    elif x7[x, 0] >= 1.5 and x7[x, 0] < 2.5:
    #elif wind67_curdifperp[x] >= 1.5 and wind67_curdifperp[x] < 2.5:
        bin_15_25_x_67 += 1
        
    else:
        bin_25_35_x_67 += 1
        
# divide by bin width and total num of points

total = [bin_n35_n25_x_67, bin_n35_n25_x_67, bin_n25_n15_x_67, bin_n25_n15_x_67,
         bin_n15_n05_x_67, bin_n15_n05_x_67, bin_n05_00_x_67, bin_00_05_x_67, 
         bin_05_15_x_67, bin_05_15_x_67, bin_15_25_x_67, bin_15_25_x_67, 
         bin_25_35_x_67, bin_25_35_x_67]



total_sum = np.sum(total)

total = total / total_sum

bins = [-3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

tlabs = ['', '', '-2.5', '', '-1.5', '', '0.5', '0.0', '0.5', '', '1.5', '',
         '2.5', '']

figure(figsize = (10, 5), dpi = 300)
plt.bar(bins, total, width = 0.5, align = 'edge', edgecolor = 'white', tick_label = tlabs)

plt.show()



