import MFT
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import cartopy as cart
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from progress_bar import update_progress

ds = Dataset('wrf_output.nc', 'r')

latdim = ds.dimensions['lat']
londim = ds.dimensions['lon']
nlat = len(latdim)
nlon = len(londim)

lat = ds['lat']
lon = ds['lon']
u10 = ds['U10']
v10 = ds['V10']
uc = ds['UC']
vc = ds['VC']
slp = ds['SLP']
t2 = ds['T2']
q2 = ds['Q2']
sst = ds['SST']

u10mask = ma.getmaskarray(u10[:])
u10_masked_indices = ma.where(u10mask == True)
u10_unmasked_indices = ma.where(u10mask == False)
ucmask = ma.getmaskarray(uc[:])
uc_masked_indices = ma.where(ucmask == True)
uc_unmasked_indices = ma.where(ucmask == False)
slpmask = ma.getmaskarray(slp[:])
slp_masked_indices = ma.where(slpmask == True)
slp_unmasked_indices = ma.where(slpmask == False)
t2mask = ma.getmaskarray(t2[:])
t2_masked_indices = ma.where(t2mask == True)
t2_unmasked_indices = ma.where(t2mask == False)
q2mask = ma.getmaskarray(q2[:])
q2_masked_indices = ma.where(q2mask == True)
q2_unmasked_indices = ma.where(q2mask == False)
sstmask = ma.getmaskarray(sst[:])
sst_masked_indices = ma.where(sstmask == True)
sst_unmasked_indices = ma.where(sstmask == False)

two_array_wind = ma.array(np.full(shape=(nlat, nlon), fill_value=2.0, dtype=float))
two_array_wind[u10_masked_indices] = ma.masked
two_array_current = ma.array(np.full(shape=(nlat, nlon), fill_value=2.0, dtype=float))
two_array_current[uc_masked_indices] = ma.masked

print('Calculating windspeeds...')

windspeed = ma.sqrt(ma.add(ma.power(u10,two_array_wind), ma.power(v10,two_array_wind)))
u_norm = ma.divide(u10,windspeed)
v_norm = ma.divide(v10,windspeed)

u_diff = ma.subtract(u10,uc)
v_diff = ma.subtract(v10,vc)
windspeed_c = ma.sqrt(ma.add(ma.power(u_diff,two_array_wind), ma.power(v_diff,two_array_wind)))
u_norm_c = ma.divide(u_diff,windspeed_c)
v_norm_c = ma.divide(v_diff,windspeed_c)

print('Calculating current speeds...')

currentspeed = ma.sqrt(ma.add(ma.power(uc,two_array_current), ma.power(vc,two_array_current)))
#uc_norm = ma.divide(uc,currentspeed)
#vc_norm = ma.divide(vc,currentspeed)

CONV_CRIT = 0.00005     #convergence critereon (fractional change)  []  
CONVECT = 0.0          #convective parameter  
warn = 1                #warning are given     
eqv_neut = 0            #output winds are winds rather than equivalent neutral winds  
z_wanted = 10.0         #height to which winds, potential temp, and humidity are adjusted                                
flux_model = 1          #1988 Smith Model
Qnet = 5.0
sst_prm = 0
z0_mom_prm = 0
z0_TQ_prm = 0
stable_prm = 0
wind_ang = 0
wave_ang = 0

dyn_in_prm = 0 # Wind speed, relative to the surface current, m/s
dyn_in_val2 = 0.0
ref_ht_wind = 10.0 # Height of the wind observations, m
CONVECT = 0 # Convective parameter. Recommended value between 0.7 and 1.25. For details see TOGA NOTES #4 (recommendation comes from no capillary waves)
air_moist_prm = 0 # Specific humidity
sfc_moist_prm = 1 # Relative humidity, fraction
sfc_moist_val = 0.98 # 98% assumed because of salinity obstructing evaporation
salinity = 34.9 / 1000.0 # Salinity, fraction, No salinity in dataset (global average)
ss_prm = 0 # Sea state parameterization, wind-wave stability parameter
ss_val = 1.0 # No sea state data in dataset, set to 1.0 for local equilibrium
ref_ht_tq = 2.0 # Height of temperature and humidity observations, m
sst_prm = 0 # Designates surface temperature as skin temperature?
astab = 1 # Atmospheric stability is calculated

water_density = 1023.6 # density of seawater in kg/m^3 at 25C, 35 g/kg salinity, and 1 atm

stress = np.ma.empty(shape = (nlat, nlon))
stress_c = np.ma.empty(shape = (nlat, nlon))

# find the indices of the dataset where ANY of the stress related variables are masked
'''
zipped_masked_indices = set(zip(u10_masked_indices[0], u10_masked_indices[1]))\
	.union(set(zip(uc_masked_indices[0], uc_masked_indices[1])))\
		.union(set(zip(slp_masked_indices[0], slp_masked_indices[1])))\
			.union(set(zip(t2_masked_indices[0], t2_masked_indices[1])))\
				.union(set(zip(q2_masked_indices[0], q2_masked_indices[1])))\
					.union(set(zip(sst_masked_indices[0], sst_masked_indices[1])))
masked_indices = list(zip(*zipped_masked_indices))
'''
masked_indices = sst_masked_indices

stress[masked_indices] = ma.masked
stress_c[masked_indices] = ma.masked
stressmask = ma.getmaskarray(stress)
unmasked_indices = np.where(stressmask == False)
unmasked_indices_count = len(unmasked_indices[0])

print('Calculating stress...')

i = 1
for (x,y) in zip(unmasked_indices[0], unmasked_indices[1]):
	update_progress(i/unmasked_indices_count)
	i += 1

	index = (x,y)
	pres = slp[index]
	pressure = pres * 100.0 # Atmospheric surface pressure, Pa (converted from hPa in dataset)
	dyn_in_val = float(windspeed[index])
	air_moist_val = float(q2[index])
	t_air = float(t2[index]) - 273.15 # Air temperature at the reference height of the thermometer and humidity sensor, C
	t_skin = float(sst[index]) - 273.15 # skin tempearture, C
	
	pass_by_ref_params = {'shf':0, 'lhf':0, 'tau':[0,0], 'u_star':[0,0], 't_star':0, 'q_star':0, 'z_over_L':0, 'wave_age':0, 'dom_phs_spd':0,
		'h_sig':0, 'ww_stab':0, 'zo_m':[0, 0], 'u_at_z':0, 't_at_z':0, 'q_at_z':0}
	count = 0
	count = MFT.ht_adj_( dyn_in_prm, dyn_in_val, dyn_in_val2, CONVECT, CONV_CRIT,
		pressure, air_moist_prm, air_moist_val, sfc_moist_prm, sfc_moist_val,
		salinity, ss_prm, ss_val, t_air, sst_prm, t_skin, ref_ht_wind, ref_ht_tq,
		z_wanted, astab, eqv_neut, Qnet, warn, flux_model, z0_mom_prm, z0_TQ_prm, stable_prm,
		pass_by_ref_params )
	if count > 0:
		stress[index] = pass_by_ref_params['tau'][0]
	else:
		stress[index] = ma.masked

	dyn_in_val = float(windspeed_c[index])
	pass_by_ref_params = {'shf':0, 'lhf':0, 'tau':[0,0], 'u_star':[0,0], 't_star':0, 'q_star':0, 'z_over_L':0, 'wave_age':0, 'dom_phs_spd':0,
		'h_sig':0, 'ww_stab':0, 'zo_m':[0, 0], 'u_at_z':0, 't_at_z':0, 'q_at_z':0}
	count = 0
	count = MFT.ht_adj_( dyn_in_prm, dyn_in_val, dyn_in_val2, CONVECT, CONV_CRIT,
		pressure, air_moist_prm, air_moist_val, sfc_moist_prm, sfc_moist_val,
		salinity, ss_prm, ss_val, t_air, sst_prm, t_skin, ref_ht_wind, ref_ht_tq,
		z_wanted, astab, eqv_neut, Qnet, warn, flux_model, z0_mom_prm, z0_TQ_prm, stable_prm,
		pass_by_ref_params )
	if count > 0:
		stress_c[index] = pass_by_ref_params['tau'][0]
	else:
		stress_c[index] = ma.masked

print('Calculating components of stress...')

ustress = ma.multiply(stress, u_norm)
vstress = ma.multiply(stress, v_norm)
ustress[masked_indices] = ma.masked
vstress[masked_indices] = ma.masked

ustress_c = ma.multiply(stress_c, u_norm_c)
vstress_c = ma.multiply(stress_c, v_norm_c)
ustress_c[masked_indices] = ma.masked
vstress_c[masked_indices] = ma.masked

print('Calculating curl of stress, Ekman net mass transport, and Ekman upwelling...')

curl = ma.masked_all(shape = (nlat, nlon), dtype=float)
curl_c = ma.masked_all(shape = (nlat, nlon))
curl_diff = ma.masked_all(shape = (nlat, nlon))

uekmantransport = ma.masked_all(shape = (nlat, nlon))
vekmantransport = ma.masked_all(shape = (nlat, nlon))
uekmantransport_c = ma.masked_all(shape = (nlat, nlon))
vekmantransport_c = ma.masked_all(shape = (nlat, nlon))
ekmanupwelling = ma.masked_all(shape = (nlat, nlon))
ekmanupwelling_c = ma.masked_all(shape = (nlat, nlon))
ekmanupwelling_diff = ma.masked_all(shape = (nlat, nlon))

omega = 7.2921 * 10 ** -5
def calc_f(latitude):
	return 2 * omega * math.sin(latitude * math.pi / 180)

i = 1
for (x,y) in zip(unmasked_indices[0], unmasked_indices[1]):
	update_progress(i/unmasked_indices_count)
	i += 1
	
	index = (x,y)

	if x != 0 and x != nlat - 1 and y != 0 and y != nlon - 1 and ustress[x,y+1] is not ma.masked and ustress[x,y-1] is not ma.masked and vstress[x+1,y] is not ma.masked and vstress[x-1,y] is not ma.masked:
		dx = (111100 * math.cos(lat[x] * math.pi / 180)) * (lon[y+1] - lon[y-1])
		dy = 111100 * (lat[x+1] - lat[x-1])
		du = ustress[x,y+1] - ustress[x,y-1]
		dv = vstress[x+1,y] - vstress[x-1,y]

		curl[index] = dv/dx - du/dy
		
		du_c = ustress_c[x,y+1] - ustress_c[x,y-1]
		dv_c = vstress_c[x+1,y] - vstress_c[x-1,y]

		curl_c[index] = dv_c/dx - du_c/dy

		curl_diff[index] = curl[index] - curl_c[index]

		f = calc_f(lat[x])

		uekmantransport[index] = vstress[index]/f
		vekmantransport[index] = -1 * ustress[index]/f
		uekmantransport_c[index] = vstress_c[index]/f
		vekmantransport_c[index] = -1 * ustress_c[index]/f

		ekmanupwelling[index] = curl[index] / (water_density * f)
		ekmanupwelling_c[index] = curl_c[index] / (water_density * f)
		ekmanupwelling_diff[index] = ekmanupwelling[index] - ekmanupwelling_c[index]

uekmantransport_diff = ma.subtract(uekmantransport, uekmantransport_c)
vekmantransport_diff = ma.subtract(vekmantransport, vekmantransport_c)
uekmantransport_diff[masked_indices] = ma.masked
vekmantransport_diff[masked_indices] = ma.masked

two_array_stress = ma.array(np.full(shape=(nlat, nlon), fill_value=2.0, dtype=float))
two_array_stress[masked_indices] = ma.masked
ekmantransport_c = ma.sqrt(ma.add(ma.power(uekmantransport_c,two_array_stress), ma.power(vekmantransport_c,two_array_stress)))
ekmantransport_diff = ma.sqrt(ma.add(ma.power(uekmantransport_diff,two_array_stress), ma.power(vekmantransport_diff,two_array_stress)))

#figsize = (10.80, 7.20)
figsize = (19.20, 10.80)
#plt.style.use('dark_background')

plot_knots = True
speed_units = r'$m s^{-1}$'
speed_scale = 50
wind_step = 4
wind_vmax = 15
current_step = 3
stress_step = 3
stress_limit = max(ma.max(stress), ma.max(stress_c))
stress_vmax = stress_limit
curl_limit = max(ma.max(ma.abs(curl)), ma.max(ma.abs(curl_c)))
curl_vmin = -1 * curl_limit
curl_vmax = curl_limit
ekmantransport_limit = ma.max(ma.abs(ekmantransport_c))
ekmantransport_vmin = 0
ekmantransport_vmax = ekmantransport_limit
ekmantransport_scale = ekmantransport_limit * 15
ekmanupwelling_limit = ma.max(ma.abs(ekmanupwelling_c))
ekmanupwelling_vmin = -1 * ekmanupwelling_limit
ekmanupwelling_vmax = ekmanupwelling_limit
extent = [float(lon[0]), float(lon[nlon-1]), float(lat[0]), float(lat[nlat-1])]

if plot_knots:
	speed_units = 'knots'
	conversion_array = np.full_like(windspeed, fill_value=1.92)
	windspeed = np.multiply(windspeed, conversion_array)
	u10 = np.multiply(u10, conversion_array)
	v10 = np.multiply(v10, conversion_array)
	currentspeed = np.multiply(currentspeed, conversion_array)
	uc = np.multiply(uc, conversion_array)
	vc = np.multiply(vc, conversion_array)
	wind_vmax *= 1.92
	speed_scale *= 1.92


fig1, axes1 = plt.subplots(nrows=2, ncols=3, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=figsize)

for index, ax in np.ndenumerate(axes1):
	ax.set_extent(extent, crs=ccrs.PlateCarree())
	ax.coastlines()
	if index != (0,0):
		ax.add_feature(cart.feature.LAND, zorder=100, edgecolor='k', facecolor='#aaa')
	gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, linestyle='--')
	gl.top_labels = False
	gl.right_labels = False
	gl.xformatter = LONGITUDE_FORMATTER
	gl.yformatter = LATITUDE_FORMATTER

axes1[0,0].set_title(f'10m Wind')
wind_contour = axes1[0,0].contourf(lon[:], lat[:], windspeed[:], 20, cmap='jet', vmin=0, vmax=wind_vmax)
axes1[0,0].figure.colorbar(wind_contour, ax=axes1[0,0]).set_label(f'wind speed ({speed_units})', rotation=-90, va='bottom')
if plot_knots:
	axes1[0,0].barbs(lon[::wind_step], lat[::wind_step], u10[::wind_step,::wind_step], v10[::wind_step,::wind_step], length=5, sizes=dict(emptybarb=0.1, spacing=0.2, height=0.5), linewidth = 0.95, transform=ccrs.PlateCarree())
else:
	axes1[0,0].quiver(lon[::wind_step], lat[::wind_step], u10[::wind_step,::wind_step], v10[::wind_step,::wind_step], scale=speed_scale*3, pivot='middle', transform=ccrs.PlateCarree())

axes1[1,0].set_title(f'Ocean Current')
current_contour = axes1[1,0].contourf(lon[:], lat[:], currentspeed[:], 20, cmap='coolwarm')
axes1[1,0].figure.colorbar(current_contour, ax=axes1[1,0]).set_label(f'current speed  ({speed_units})', rotation=-90, va='bottom')
axes1[1,0].quiver(lon[::current_step], lat[::current_step], uc[::current_step,::current_step], vc[::current_step,::current_step], scale=speed_scale, pivot='middle', transform=ccrs.PlateCarree())

axes1[0,1].set_title(f'Wind Stress (ignoring currents)')
stress_contour = axes1[0,1].contourf(lon[:], lat[:], stress[:], 20, cmap='plasma', vmin=0, vmax=stress_vmax)
axes1[0,1].figure.colorbar(stress_contour, ax=axes1[0,1]).set_label(f'wind stress  ({r"$N m^{-2}$"})', rotation=-90, va='bottom')
axes1[0,1].quiver(lon[::stress_step], lat[::stress_step], ustress[::stress_step,::stress_step], vstress[::stress_step,::stress_step], pivot='middle', transform=ccrs.PlateCarree())

axes1[1,1].set_title(f'Wind Stress')
stress_c_contour = axes1[1,1].contourf(lon[:], lat[:], stress_c[:], 20, cmap='plasma', vmin=0, vmax=stress_vmax)
axes1[1,1].figure.colorbar(stress_c_contour, ax=axes1[1,1]).set_label(f'wind stress  ({r"$N m^{-2}$"})', rotation=-90, va='bottom')
axes1[1,1].quiver(lon[::stress_step], lat[::stress_step], ustress_c[::stress_step,::stress_step], vstress_c[::stress_step,::stress_step], pivot='middle', transform=ccrs.PlateCarree())

axes1[0,2].set_title(f'Vectors: Wind Stress (ignoring currents), Shading: Curl')
stresscurl_contour = axes1[0,2].contourf(lon[:], lat[:], curl[:], 20, cmap='bwr', vmin=curl_vmin, vmax=curl_vmax)
axes1[0,2].figure.colorbar(stresscurl_contour, ax=axes1[0,2]).set_label(f'curl of wind stress ({r"$N m^{-2}m^{-1}$"})', rotation=-90, va='bottom')
axes1[0,2].quiver(lon[::stress_step], lat[::stress_step], ustress[::stress_step,::stress_step], vstress[::stress_step,::stress_step], pivot='middle', transform=ccrs.PlateCarree())

axes1[1,2].set_title(f'Vectors: Wind Stress, Shading: Curl of Wind Stress')
stresscurl_c_contour = axes1[1,2].contourf(lon[:], lat[:], curl_c[:], 20, cmap='bwr', vmin=curl_vmin, vmax=curl_vmax)
axes1[1,2].figure.colorbar(stresscurl_c_contour, ax=axes1[1,2]).set_label(f'curl of wind stress ({r"$N m^{-2}m^{-1}$"})', rotation=-90, va='bottom')
axes1[1,2].quiver(lon[::stress_step], lat[::stress_step], ustress_c[::stress_step,::stress_step], vstress_c[::stress_step,::stress_step], pivot='middle', transform=ccrs.PlateCarree())


fig2, axes2 = plt.subplots(nrows=2, ncols=3, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=figsize)

for index, ax in np.ndenumerate(axes2):
	ax.set_extent(extent, crs=ccrs.PlateCarree())
	ax.coastlines()
	ax.add_feature(cart.feature.LAND, zorder=100, edgecolor='k', facecolor='#aaa')
	gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, linestyle='--')
	gl.top_labels = False
	gl.right_labels = False
	gl.xformatter = LONGITUDE_FORMATTER
	gl.yformatter = LATITUDE_FORMATTER

axes2[0,0].set_title(f'Vectors: Wind Stress, Shading: Curl of Wind Stress')
stresscurl_c_contour2 = axes2[0,0].contourf(lon[:], lat[:], curl_c[:], 20, cmap='bwr', vmin=curl_vmin, vmax=curl_vmax)
axes2[0,0].figure.colorbar(stresscurl_c_contour2, ax=axes2[0,0]).set_label(f'curl of wind stress ({r"$N m^{-2}m^{-1}$"})', rotation=-90, va='bottom')
axes2[0,0].quiver(lon[::stress_step], lat[::stress_step], ustress_c[::stress_step,::stress_step], vstress_c[::stress_step,::stress_step], pivot='middle', transform=ccrs.PlateCarree())

axes2[1,0].set_title(f'Error in Wind Stress Curl')
stresscurl_diff_contour = axes2[1,0].contourf(lon[:], lat[:], curl_diff[:], 20, cmap='bwr', vmin=curl_vmin, vmax=curl_vmax)
axes2[1,0].figure.colorbar(stresscurl_diff_contour, ax=axes2[1,0]).set_label(f'curl of wind stress ({r"$N m^{-2}m^{-1}$"})', rotation=-90, va='bottom')

axes2[0,1].set_title(f'Ekman Net Mass Transport')
ekmantransport_c_contour = axes2[0,1].contourf(lon[:], lat[:], ekmantransport_c[:], 20, cmap='plasma', vmin=ekmantransport_vmin, vmax=ekmantransport_vmax)
axes2[0,1].figure.colorbar(ekmantransport_c_contour, ax=axes2[0,1]).set_label(f'mass transport ({r"$kg s^{-1}m^{-1}$"})', rotation=-90, va='bottom')
axes2[0,1].quiver(lon[::stress_step], lat[::stress_step], uekmantransport_c[::stress_step,::stress_step], vekmantransport_c[::stress_step,::stress_step], scale=ekmantransport_scale, pivot='middle', transform=ccrs.PlateCarree())

axes2[1,1].set_title(f'Error in Ekman Net Mass Transport')
ekmantransport_diff_contour = axes2[1,1].contourf(lon[:], lat[:], ekmantransport_diff[:], 20, cmap='plasma', vmin=ekmantransport_vmin, vmax=ekmantransport_vmax)
axes2[1,1].figure.colorbar(ekmantransport_diff_contour, ax=axes2[1,1]).set_label(f'mass transport ({r"$kg s^{-1}m^{-1}$"})', rotation=-90, va='bottom')
axes2[1,1].quiver(lon[::stress_step], lat[::stress_step], uekmantransport_diff[::stress_step,::stress_step], vekmantransport_diff[::stress_step,::stress_step], scale=ekmantransport_scale, pivot='middle', transform=ccrs.PlateCarree())

axes2[0,2].set_title(f'Ekman Upwelling')
ekmanupwelling_c_contour = axes2[0,2].contourf(lon[:], lat[:], ekmanupwelling_c[:], 20, cmap='bwr', vmin=ekmanupwelling_vmin, vmax=ekmanupwelling_vmax)
axes2[0,2].figure.colorbar(ekmanupwelling_c_contour, ax=axes2[0,2]).set_label(f'vertical velocity ({r"$ms^{-1}$"})', rotation=-90, va='bottom')

axes2[1,2].set_title(f'Error in Ekman Upwelling')
ekmanupwelling_diff_contour = axes2[1,2].contourf(lon[:], lat[:], ekmanupwelling_diff[:], 20, cmap='bwr', vmin=ekmanupwelling_vmin, vmax=ekmanupwelling_vmax)
axes2[1,2].figure.colorbar(ekmanupwelling_diff_contour, ax=axes2[1,2]).set_label(f'vertical velocity ({r"$ms^{-1}$"})', rotation=-90, va='bottom')

t2_conversion_array = np.full_like(t2, fill_value=273.15)
sst_conversion_array = np.full_like(sst, fill_value=273.15)

fig3, axes3 = plt.subplots(nrows=2, ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=figsize)

for index, ax in np.ndenumerate(axes3):
	ax.set_extent(extent, crs=ccrs.PlateCarree())
	ax.coastlines()
	if index != (0,0) and index != (0,0):
		ax.add_feature(cart.feature.LAND, zorder=100, edgecolor='k', facecolor='#aaa')
	gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, linestyle='--')
	gl.top_labels = False
	gl.right_labels = False
	gl.xformatter = LONGITUDE_FORMATTER
	gl.yformatter = LATITUDE_FORMATTER

axes3[0,0].set_title(f'10m Wind (barbs and shading) and SLP (contours)')
wind_contour = axes3[0,0].contourf(lon[:], lat[:], windspeed[:], 20, cmap='jet', vmin=0, vmax=wind_vmax)
axes3[0,0].figure.colorbar(wind_contour, ax=axes3[0,0]).set_label(f'wind speed ({speed_units})', rotation=-90, va='bottom')
pressure_contour = axes3[0,0].contour(lon[:], lat[:], slp[:], colors=['white'], levels=5)
axes3[0,0].clabel(pressure_contour, colors=['white'], manual=False, inline=True, fmt=' {:.0f} '.format)
if plot_knots:
	axes3[0,0].barbs(lon[::wind_step], lat[::wind_step], u10[::wind_step,::wind_step], v10[::wind_step,::wind_step], length=5, sizes=dict(emptybarb=0.1, spacing=0.2, height=0.5), linewidth = 0.95)
else:
	axes3[0,0].quiver(lon[::wind_step], lat[::wind_step], u10[::wind_step,::wind_step], v10[::wind_step,::wind_step], scale=speed_scale*3, pivot='middle')

air_sea_temp_diff = ma.subtract(t2, sst)
air_sea_temp_diff = ma.masked_where(air_sea_temp_diff < -3, air_sea_temp_diff)
air_sea_temp_diff = ma.masked_where(air_sea_temp_diff > 3, air_sea_temp_diff)

axes3[0,1].set_title(f'Air-Sea Temperature Difference')
t2_contour = axes3[0,1].contourf(lon[:], lat[:], air_sea_temp_diff[:], 20, cmap='bwr', vmin=-3, vmax=3)
axes3[0,1].figure.colorbar(t2_contour, ax=axes3[0,1]).set_label(f'temperature (C)', rotation=-90, va='bottom')

axes3[1,0].set_title(f'Ocean Current')
current_contour = axes3[1,0].contourf(lon[:], lat[:], currentspeed[:], 20, cmap='coolwarm')
axes3[1,0].figure.colorbar(current_contour, ax=axes3[1,0]).set_label(f'current speed ({speed_units})', rotation=-90, va='bottom')
axes3[1,0].quiver(lon[::current_step], lat[::current_step], uc[::current_step,::current_step], vc[::current_step,::current_step], scale=speed_scale, pivot='middle')

sst_celcius = ma.subtract(sst, sst_conversion_array)
sst_celcius = ma.masked_where(sst_celcius < 15, sst_celcius)

axes3[1,1].set_title(f'Sea Surface Temperature')
sst_contour = axes3[1,1].contourf(lon[:], lat[:], sst_celcius[:], 20, cmap='coolwarm', vmin=15, vmax=30)
axes3[1,1].figure.colorbar(sst_contour, ax=axes3[1,1]).set_label(f'temperature (C)', rotation=-90, va='bottom')
axes3[1,1].quiver(lon[::current_step], lat[::current_step], uc[::current_step,::current_step], vc[::current_step,::current_step], scale=speed_scale, pivot='middle')


plt.show()
ds.close()