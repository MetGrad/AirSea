import math
import xarray as xr
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt 
import matplotlib.pylab as plt
import matplotlib.path as mpath
import numpy as np
import cartopy.feature as cfeature
import metpy.calc as mpcalc
from metpy.units import units
import pandas as pd 

#Varibales below on 0.08 degree grid in Gulf of Mexico 

# 2 metre dewpoint temperature in K 
#(initial_time0_hours: 4, g4_lat_1: 185, g4_lon_2: 299)
#2D_GDS4_SFC
dsTdew = xr.open_dataset('/Users/Anna/Desktop/MSMET/AirSea/Project1/ECMWFsfc/ec.oper.an.sfc.128_168_2d.regn1280sc.20180326.grb.smoot576067.nc')
#print(dsTdew.keys)
#print(dsTdew.variables)

#2 metre temperature in K 
#(initial_time0_hours: 4, g4_lat_1: 185, g4_lon_2: 299)
#2T_GDS4_SFC
dsT = xr.open_dataset('/Users/Anna/Desktop/MSMET/AirSea/Project1/ECMWFsfc/ec.oper.an.sfc.128_167_2t.regn1280sc.20180326.grb.smoot576067.nc')
#print(dsT.keys)
#print(dsT.variables)

#U
#10 metre U wind component in m/2
#(initial_time0_hours: 4, g4_lat_1: 185, g4_lon_2: 299)
#'10U_GDS4_SFC'
dsU = xr.open_dataset('/Users/Anna/Desktop/MSMET/AirSea/Project1/ECMWFsfc/ec.oper.an.sfc.128_165_10u.regn1280sc.20180326.grb.smoot576067.nc')
#print(dsU.keys)
#print(dsU.variables)

#V
#10 metre V wind component in m/s
#(initial_time0_hours: 4, g4_lat_1: 185, g4_lon_2: 299
#10V_GDS4_SFC
dsV = xr.open_dataset('/Users/Anna/Desktop/MSMET/AirSea/Project1/ECMWFsfc/ec.oper.an.sfc.128_166_10v.regn1280sc.20180326.grb.smoot576067.nc')
#print(dsV.keys)
#print(dsV.variables)

#Surface pressure in Pa
##(initial_time0_hours: 4, g4_lat_1: 185, g4_lon_2: 299
#SP_GDS4_SFC
dsP = xr.open_dataset('/Users/Anna/Desktop/MSMET/AirSea/Project1/ECMWFsfc/ec.oper.an.sfc.128_134_sp.regn1280sc.20180326.grb.smoot576067.nc')
#print(dsP.keys)
#print(dsP.variables)

#Skin Temp is SST

#'SKT_GDS4_SFC'
#g4_lon_2
#Dimensions:(initial_time0_hours: 4, g4_lat_1: 185, g4_lon_2: 299)
#Coordinates:
#  * g4_lat_1               (g4_lat_1) float32 30.97 30.9 30.83 ... 18.1 18.03
#  * g4_lon_2               (g4_lon_2) float32 261.0 261.1 261.1 ... 281.9 282.0
#  * initial_time0_hours    (initial_time0_hours) datetime64[ns] 2018-03-26 .....
# Time: (['2018-03-26T00:00:00.000000000', '2018-03-26T06:00:00.000000000',
       #'2018-03-26T12:00:00.000000000', '2018-03-26T18:00:00.000000000']

dsSkin = xr.open_dataset('/Users/Anna/Desktop/MSMET/AirSea/Project1/ECMWFsfc/ec.oper.an.sfc.128_235_skt.regn1280sc.20180326.grb.smoot576067.nc')

#print(ds.keys)
#print(ds.variables)

SkinT=dsSkin.SKT_GDS4_SFC.values[0,:,:] #in K @ sm
SkinTC=SkinT-273.15 #in C 

AtmosT = dsT['2T_GDS4_SFC'].values[0,:,:]
#AtmosT=dsT.2T_GDS4_SFC.values[0,:,:] #in K @ 2m
AtmosTC=AtmosT-273.15 #in C 

#has zeros 
Uwind=dsU['10U_GDS4_SFC'].values[0,:,:] #in m/s @ 10m
Vwind= dsV['10V_GDS4_SFC'].values[0,:,:] #in m/s @ 10m
windMag = np.sqrt(Uwind**2+Vwind**2) 

sfcP= dsP.SP_GDS4_SFC.values[0,:,:] #in Pa 

#May have zeros 
Tdew=dsTdew['2D_GDS4_SFC'].values[0,:,:] #in K @ 2m
TdewC=Tdew-273.15 #in C 
   
Lat=dsSkin.SKT_GDS4_SFC.g4_lat_1.values
Lon=dsSkin.SKT_GDS4_SFC.g4_lon_2.values

#(initial_time0_hours: 4, g4_lat_1: 185, g4_lon_2: 299) 
xLen=len(Lon)
yLen=len(Lat)

print(Lat.shape)
print(Lon.shape)

#PressurePlot 

mapcrs = ccrs.PlateCarree()

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection=mapcrs)

#gives coastlines
ax.coastlines()
#filled contours (lon, lat, thing plotted)
cs = ax.contourf(Lon, Lat, sfcP)
#contnour lines
ax.contour(Lon, Lat, sfcP, colors = 'white', linewidths = 1, transform=ccrs.PlateCarree())
#ax=plt.gca()
PCM=ax.get_children()[2]
plt.colorbar(cs, ax=ax, shrink = 0.5, label = 'Pa')
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, color = 'black')
#don't want labels
gl.top_labels = False
gl.ylabels_right = False
plt.title('Pressure')
#plt.savefig('PV_110.png', bbox_inches='tight', dpi=150)
plt.show()

from __future__ import print_function
from MFT21 import *
import sys

CONV_CRIT = 0.00005     #convergence critereon (fractional change)  []  
CONVECT = 0.0           #convective parameter  
warn = 0                #warning are given for 1, and hidden for 0    
eqv_neut_prm = 0        #output winds are winds rather than equivalent neutral winds  
z_wanted = 10.0         #height to which winds, potential temp, and humidity are adjusted                                
flux_model = -1         #BVW model=0  

Qnet = 5.0
sst_prm = 2.0 #wave age
z0_mom_prm = 0 #like the smith 88 model 
z0_TQ_prm = 1 #0, 1, 4, 5 are all ok 
#Stability Parameter 0=BusingerDyer,1=SmithCollection,2=COAWST,3=Paulson1970
#3 was default
stable_prm = 3
A_oil = 0.0

#Inputs I set 
dyn_in_prm = 0.0
dyn_in_val2 = 0.0
air_moist_prm = 2
ref_ht_wind = 10.0
ref_ht_tq = 2.0
sfc_moist_prm = 1.0
sfc_moist_val = 0.98 #over water 
salinity = 0.349 #global average
ss_prm = 0.0 
ss_val = 28 #42.8 #typical of ocean/wind swell
astab = 1.0

#Outputs 
#shf = 0.0
#lhf = 0.0
#tau = [0.0,0.0]
#u_star = [0.0,0.0]
#t_star = 0.0
#q_star= 0.0
#z_over_L = 0.0
#wave_age = 0.0
#dom_phase_spd = 0.0
#hsig = 0.0
#ww_stab = 0.0 
#zo_m = 0.0001
#u_at_z = 0.0
#t_at_z = 0.0 
#q_at_z = 0.0



SkinT=dsSkin.SKT_GDS4_SFC.values[0,:,:] #in K @ sm
SkinTC=SkinT-273.15 #in C 

AtmosT = dsT['2T_GDS4_SFC'].values[0,:,:]
AtmosTC=AtmosT-273.15 #in C 

Uwind=dsU['10U_GDS4_SFC'].values[0,:,:] #in m/s @ 10m
Vwind= dsV['10V_GDS4_SFC'].values[0,:,:] #in m/s @ 10m
windMag = np.sqrt(Uwind**2+Vwind**2) 

sfcP= dsP.SP_GDS4_SFC.values[0,:,:] #in Pa 

Tdew=dsTdew['2D_GDS4_SFC'].values[0,:,:] #in K @ 2m
TdewC=Tdew-273.15 #in C 
   
Lat=dsSkin.SKT_GDS4_SFC.g4_lat_1.values
Lon=dsSkin.SKT_GDS4_SFC.g4_lon_2.values

xLen=len(Lon)
yLen=len(Lat)


#MFT Data

#fill new arrays with zeros here 
shfArr3=np.zeros((yLen,xLen))
lhfArr3=np.zeros((yLen,xLen))
tauArr3=np.zeros((yLen,xLen))
uStarArr3=np.zeros((yLen,xLen))
tStarArr3=np.zeros((yLen,xLen))
qStarArr3=np.zeros((yLen,xLen))
zOverlArr3=np.zeros((yLen,xLen))
waveAgeArr3=np.zeros((yLen,xLen))
domPhaseSpeedArr3=np.zeros((yLen,xLen))
hsigArr3=np.zeros((yLen,xLen))
wwStabArr3=np.zeros((yLen,xLen))
zomArr3=np.zeros((yLen,xLen))
uAtzArr3=np.zeros((yLen,xLen))
tAtzArr3=np.zeros((yLen,xLen))
qAtzArr3=np.zeros((yLen,xLen))

#fill each array with nans 
shfArr3=np.nan
lhfArr3=np.nan
tauArr3=np.nan
uStarArr3=np.nan
tStarArr3=np.nan
qStarArr3=np.nan
zOverlArr3=np.nan
waveAgeArr3=np.nan
domPhaseSpeedArr3=np.nan
hsigArr3=np.nan
wwStabArr3=np.nan
zomArr3=np.nan
uAtzArr3=np.nan
tAtzArr3=np.nan
qAtzArr3=np.nan


for x in range(0,xLen):
    for y in range(0,yLen):
        #assign other input arrays here
        dyn_in_val = windMag[y,x]
        pressure = sfcP[y,x]
        air_moist_val = TdewC[y,x]
        t_air=AtmosTC[y,x]
        t_skin=SkinTC[y,x]
        
        if dyn_in_val==0.0:
            continue 
        if pressure==0.0:
            continue
        if air_moist_val==0.0:
            continue
        if t_air==0.0:
            continue
        if t_skin==0.0:
            continue
                
        try:

            count, shf, lhf, tau, u_star, t_star, q_star, z_over_L, wave_age, dom_phase_spd, hsig, ww_stab, zo_m, u_at_z, t_at_z, q_at_z = ht_adj_( dyn_in_prm, 
                    dyn_in_val, dyn_in_val2, CONVECT, CONV_CRIT,
                    pressure, air_moist_prm, air_moist_val, sfc_moist_prm, sfc_moist_val,
                    salinity, ss_prm, ss_val, t_air, sst_prm, t_skin, ref_ht_wind, ref_ht_tq,
                    z_wanted, astab, eqv_neut_prm, Qnet, warn, flux_model, z0_mom_prm, z0_TQ_prm, stable_prm,
                    A_oil, u_star, t_star, q_star, z_over_L, wave_age, dom_phase_spd, hsig, ww_stab, zo_m, u_at_z, t_at_z, q_at_z )

        except ZeroDivisionError:
            #print whatever type of error it is, and print the values associated with the error 
            print("Problem!")
            print(dyn_in_val)
            print(pressure)
            print(air_moist_val)
            print(t_air)
            print(t_skin)
            

        #shf as output from function above
        #assign other array outputs here
        shfArr3[y,x]=shf
        lhfArr3[y,x]=lhf
        tauArr3[y,x]=tau
        uStar3[y,x]=u_star
        tStarArr3[y,x]=t_star
        qStarArr3[y,x]=q_star
        zOverlArr3[y,x]=z_over_l
        waveAgeArr3[y,x]=wave_age
        domPhaseSpeedArr3[y,x]=dom_phase_spd
        hsigArr3[y,x]=hsig
        wwStabArr3[y,x]=ww_stab
        zomArr3[y,x]=zo_m
        uAtzArr3[y,x]=u_at_z
        tAtzArr3[y,x]=t_at_z
        qAtzArr3[y,x]=q_at_z


