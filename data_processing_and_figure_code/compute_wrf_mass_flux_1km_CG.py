#==========================================================================
# Title: compute_wrf_mass_flux_1km_CG.py
# Author: McKenna W. Stanford
# Utility: Reads in coarse-grained 1-km simulation and computes mass flux & TWC.
# Writes these variables along with vertical velocity (w)
# to a dictionary named 'baseline_1km_CG_mass_flux_twc_w_dict.p'
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
import pandas as pd



path = '/glade/scratch/mckenna/AMIE/baseline_1km/3km/'

files = sorted(glob.glob(path+'/wrfout*'))
nt = len(files)
    

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
    
    ncfile = xarray.open_dataset(files[tt])
    w = ncfile['w'].values
    rho_air = ncfile['rho_air'].values
    twc = ncfile['twc'].values
    lon = ncfile['lon'].values
    lat = ncfile['lat'].values
    z = ncfile['z'].values
    time = ncfile['time'].values
    time = pd.to_datetime(time,unit='s')
    ncfile.close()

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
time_arr = time_arr[:,0]

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
f = open(savdir+'baseline_1km_CG_mass_flux_twc_w_dict.p','wb')
pickle.dump(out_dict,f)
f.close()

print('Completed writing to dictionary. Done.')


