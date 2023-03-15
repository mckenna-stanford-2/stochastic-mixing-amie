#==========================================================================
# Title: compute_wrf_mass_flux_dicts.py
# Author: McKenna W. Stanford
# Utility: Reads in wrfout* files and computes the mass flux and twc. 
#Also
# vertically integrates mass flux considering only upward or downward motions
# and upward or downward motions conditionally sampled for grid points with
# vertical velocity > 1 m/s and for veertical velocity > 2 m/s (i.e., convective
# updrafts with differing thresholds.)
# Writes these variables along with vertical velocity (w)
# to a dictionary named '<SIM_NAME>_mass_flux_twc_w_dict.p'
#==========================================================================

#------------------------------------------
# Imports
#------------------------------------------
import numpy as np
import glob
import xarray
import matplotlib
import pickle
import pandas as pd
from netCDF4 import Dataset
import wrf

sim_names = ['baseline',\
             'stoch1_short',\
             'stoch2_short',\
             'stoch3_short',\
             'stoch4_short',\
             'stoch5_short',\
             'stoch1_long',\
             'stoch2_long',\
             'stoch3_long',\
             'stoch4_long',\
             'stoch5_long',\
             'theta_pert1',\
             'theta_pert2',\
             'theta_pert3',\
             'theta_pert4',\
             'theta_pert5',\
             'no_mixing',\
             '4x',\
            ]

#sim_names = sim_names[0:8]
sim_names = sim_names[8:]
#sim_names = sim_names[8:12
#sim_names = sim_names[12:16]
#sim_names = sim_names[16:]

num_sims = len(sim_names)

path = '/glade/scratch/mckenna/AMIE/'
#=========================================================
# WRF vertical mass flux, twc, and w
#=========================================================
print('Computing mass flux...')
for sim_ii in range(num_sims):
    print('Simulation:',sim_names[sim_ii])
    files = sorted(glob.glob(path+sim_names[sim_ii]+'/wrfout*'))
    
    # limit to 12 hours after spinup
    files = files[48:]
    nt = len(files)
    #nt=3
    
    # Append to list as we loop through time
    mass_flux_arr = []
    w_arr = []
    z_arr = []
    twc_arr = []
    time_arr = []
    rho_air_arr = []

    #========================================================
    # Loop through times
    #========================================================
    for tt in range(nt):
        print('Time step {}/{}:'.format(tt+1,nt))  
        
        ncfile = Dataset(files[tt])
        lat = wrf.getvar(ncfile,'lat',meta=False).data.T
        lon = wrf.getvar(ncfile,'lon',meta=False).data.T
        w = wrf.getvar(ncfile,'wa',meta=False).data.T # units of m/s
        qr = wrf.getvar(ncfile,'QRAIN',meta=False).data.T # units of kg/kg 
        qc = wrf.getvar(ncfile,'QCLOUD',meta=False).data.T # units of kg/kg
        qi = wrf.getvar(ncfile,'QICE',meta=False).data.T # units of kg/kg
        tv = wrf.getvar(ncfile,'tv',units='K',meta=False).data.T # virtual potential temperature, units of K
        pres = wrf.getvar(ncfile,'pres',units='hPa',meta=False).data.T # pressure, units of hPa
        z = wrf.getvar(ncfile,'z',meta=False).data.T # height AGL, units of m
        time = pd.to_datetime(wrf.getvar(ncfile,'times').data)
        ncfile.close()
        
        qt = qi+qr+qc # convert total condensate mixing ratio
        rd = 287.04
        # compute air density
        rho_air = 100.*pres/(rd*tv)
        # Compute TWC and convert to g/m^3
        twc = rho_air*qt*1.e3
        # Zero out small negative values
        twc[twc < 0.] = 0.
        # compute mass flux as product of vertical velocity and and air density
        # NOTE! Does not include condensate flux
        mass_flux = w*rho_air
        
        mass_flux_arr.append(mass_flux)
        z_arr.append(z)
        rho_air_arr.append(rho_air)
        w_arr.append(w)
        twc_arr.append(twc)
        time_arr.append(time)
                
    mass_flux_arr = np.array(mass_flux_arr)
    w_arr = np.array(w_arr)
    z_arr = np.array(z_arr)
    twc_arr = np.array(twc_arr)
    rho_air_arr = np.array(rho_air_arr)
    time_arr = np.array(time_arr)


    out_dict = {'mass_flux':mass_flux_arr,\
                'w':w_arr,\
                'twc':twc_arr,\
                'rho_air':rho_air_arr,\
                'lon':lon,\
                'lat':lat,\
                'z':z_arr,\
                'time':time_arr}
    
    # Print diagnostics
    if False:
    #if True:
        for key,val in out_dict.items():
            if key != 'time':
                print(key,np.shape(val))
            else:
                print(key,val)
                
    print('Writing to dictionary...')

    savdir = '/glade/scratch/mckenna/AMIE/amie_post/mass_flux/'
    f = open(savdir+sim_names[sim_ii]+'_mass_flux_twc_w_dict.p','wb')
    pickle.dump(out_dict,f)
    f.close()
    
    print('Completed writing to dictionary. Done.')
    

