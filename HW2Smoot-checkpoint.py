from __future__ import print_function
from MFT import *
import sys

CONV_CRIT = 0.00005     #convergence critereon (fractional change)  []  
CONVECT = 0.0          #convective parameter  
warn = 1                #warning are given     
eqv_neut = 0            #output winds are winds rather than equivalent neutral winds  
z_wanted = 10.0         #height to which winds, potential temp, and humidity are adjusted                                
flux_model = 9          #BVW model=0  
Qnet = 5.0
sst_prm = 0
z0_mom_prm = 0
z0_TQ_prm = 0
stable_prm = 0
wind_ang = 0
wave_ang = 0


pass_by_ref_params = {'shf':0, 'lhf':0, 'tau':[0,0], 'u_star':[0,0], 't_star':0, 'q_star':0, 'z_over_L':0, 'wave_age':0, 'dom_phs_spd':0,
                        'h_sig':0, 'ww_stab':0, 'zo_m':[0, 0], 'u_at_z':0, 't_at_z':0, 'q_at_z':0}
    
try:
    #data_file = open("testdatav2.dat", "r")
    data_file = open("testdata12.dat", "r")
except IOError:
    print("Could not read file.")
    sys.exit()
        
data_file.readline()
print( "run  U   |ustar| ustar1 ustar2 tstar   qstar     zref/L     cp     wa   Hsig   tau1   tau2    shf   lhf    u(z)   t(z)   q(z)")
#change 62 to 5 (0-5)
for x in range(0,62):
    data_in_row = data_file.readline().split()
    for a in [2, 3, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16]:
        data_in_row[a] = float(data_in_row[a])
    for b in [0, 1, 4, 6, 8, 17]:
        data_in_row[b] = int(data_in_row[b])
        
    if len(data_in_row) == 18 :  
        num, dyn_in_prm, dyn_in_val, dyn_in_val2, ss_prm, ss_val, air_moist_prm, air_moist_val, sfc_moist_prm, sfc_moist_val, t_skin, t_air, ref_ht_wind, ref_ht_tq, pressure, salinity, CONVECT, astab = data_in_row
        wind_ang = wind_ang - 270.0;
        wave_ang = wave_ang - 90;
            
        count = ht_adj_( dyn_in_prm, dyn_in_val, dyn_in_val2, CONVECT, CONV_CRIT,
                pressure, air_moist_prm, air_moist_val, sfc_moist_prm, sfc_moist_val,
                salinity, ss_prm, ss_val, t_air, sst_prm, t_skin, ref_ht_wind, ref_ht_tq,
                z_wanted, astab, eqv_neut, Qnet, warn, flux_model, z0_mom_prm, z0_TQ_prm, stable_prm,
                pass_by_ref_params );
        
        if count <= 1 :
            print("non-convergence:" )

        print( "%2i %5.2f %6.3f %6.3f %6.3f %7.4f %9.6f %8.5f %6.2f %6.2f %5.2f %6.3f %6.3f %6.2f %6.2f %6.2f %6.2f %6.3f" %
            (x+1, dyn_in_val, m.sqrt( pass_by_ref_params['u_star'][0] * pass_by_ref_params['u_star'][0] + pass_by_ref_params['u_star'][1] * pass_by_ref_params['u_star'][1] ),
            pass_by_ref_params['u_star'][0], pass_by_ref_params['u_star'][1], pass_by_ref_params['t_star'], pass_by_ref_params['q_star'],
            pass_by_ref_params['z_over_L'], pass_by_ref_params['dom_phs_spd'],
            pass_by_ref_params['wave_age'], pass_by_ref_params['h_sig'], pass_by_ref_params['tau'][0], pass_by_ref_params['tau'][1],
            pass_by_ref_params['shf'], pass_by_ref_params['lhf'], pass_by_ref_params['u_at_z'], pass_by_ref_params['t_at_z'],
            pass_by_ref_params['q_at_z']), sep="")
