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
sst_prm = 0
z0_mom_prm = 7
z0_TQ_prm = 1
stable_prm = 3
A_oil = 0.0
shf = 0.0
lhf = 0.0
tau = [0.0,0.0]
u_star = [0.0,0.0]
t_star = 0.0
q_star= 0.0
z_over_L = 0.0
wave_age = 0.0
dom_phase_spd = 0.0
hsig = 0.0
ww_stab = 0.0 
zo_m = 0.0001
u_at_z = 0.0
t_at_z = 0.0 
q_at_z = 0.0

    
try:
    data_file = open("testdata_MFT.dat", "r")
except IOError:
    print("Could not read file.")
    sys.exit()
        
data_file.readline()
print( "run  U   |ustar| ustar1 ustar2 tstar   qstar     zref/L     cp     wa   Hsig   tau1   tau2    shf   lhf    u(z)   t(z)   q(z)")
for x in range(0,62):
    data_in_row = data_file.readline().split()
    for a in [2, 3, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16]:
        data_in_row[a] = float(data_in_row[a])
    for b in [0, 1, 4, 6, 8, 17]:
        data_in_row[b] = int(data_in_row[b])
        
    if len(data_in_row) == 18 :  
        num, dyn_in_prm, dyn_in_val, dyn_in_val2, ss_prm, ss_val, air_moist_prm, air_moist_val, sfc_moist_prm, sfc_moist_val, t_skin, t_air, ref_ht_wind, ref_ht_tq, pressure, salinity, CONVECT, astab = data_in_row
# two lines of added code to show now to keep z/L equal to an input value. Note 'z' in z/L is the wind reference height.    
             
        count, shf, lhf, tau, u_star, t_star, q_star, z_over_L, wave_age, dom_phase_spd, hsig, ww_stab, zo_m, u_at_z, t_at_z, q_at_z = ht_adj_( dyn_in_prm, 
                dyn_in_val, dyn_in_val2, CONVECT, CONV_CRIT,
                pressure, air_moist_prm, air_moist_val, sfc_moist_prm, sfc_moist_val,
                salinity, ss_prm, ss_val, t_air, sst_prm, t_skin, ref_ht_wind, ref_ht_tq,
                z_wanted, astab, eqv_neut_prm, Qnet, warn, flux_model, z0_mom_prm, z0_TQ_prm, stable_prm,
                A_oil, u_star, t_star, q_star, z_over_L, wave_age, dom_phase_spd, hsig, ww_stab, zo_m, u_at_z, t_at_z, q_at_z )

        #if count <= 1 :
        #    print("non-convergence:" )

        print( "%2i %5.2f %6.3f %6.3f %6.3f %7.4f %9.6f %8.5f %6.2f %6.2f %5.2f %6.3f %6.3f %6.2f %6.2f %6.2f %6.2f %6.3f %8.6f" %
            (x+1, dyn_in_val, m.sqrt( u_star[0] * u_star[0] + u_star[1] * u_star[1] ), u_star[0], u_star[1], t_star, q_star,
            z_over_L, dom_phase_spd, wave_age, hsig, tau[0], tau[1], shf, lhf, u_at_z, t_at_z, q_at_z,  zo_m), sep="")


