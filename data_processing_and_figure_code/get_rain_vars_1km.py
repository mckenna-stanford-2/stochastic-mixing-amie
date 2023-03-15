#======================================================================
# Title: get_rain_vars_1km.py
# Author: McKenna W. Stanford
# Utility: Reads in wrfout* files for the 1-km simulation and grabs 
# variables (qr,nr,and rho_air) that are needed to compute instantaneous 
# rain rates from P3.
# Writes to a dictionary named baseline_1km_rain_vars_at_2.5-km.p'
#======================================================================

import numpy as np
import pickle
import xarray
import scipy.special as scs
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature
from matplotlib.cm import get_cmap
import matplotlib
import glob
from netCDF4 import Dataset
import wrf
import pickle
import pandas as pd

#------------------------------------------
# Functions
#------------------------------------------
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

path = '/glade/scratch/mckenna/AMIE/baseline_1km/'

files = sorted(glob.glob(path+'/wrfout*'))
nt = len(files)

qr_arr= []
nr_arr = []
rho_air_arr = []
z_arr = []
time_arr = []

for tt in range(nt):
    print('Time step {}/{}:'.format(tt+1,nt))  
    ncfile = Dataset(files[tt])
    if tt == 0.:
        lat = wrf.getvar(ncfile,'lat',meta=False).data.T
        lon = wrf.getvar(ncfile,'lon',meta=False).data.T
    qr = wrf.getvar(ncfile,'QRAIN',meta=False).data.T # units of kg/kg 
    nr = wrf.getvar(ncfile,'QNRAIN',meta=False).data.T # units of /kg 
    tv = wrf.getvar(ncfile,'tv',units='K',meta=False).data.T # virtual potential temperature, units of K
    pres = wrf.getvar(ncfile,'pres',units='hPa',meta=False).data.T # pressure, units of hPa
    z = wrf.getvar(ncfile,'z',meta=False).data.T # height AGL, units of m
    time = pd.to_datetime(wrf.getvar(ncfile,'times').data)
    ncfile.close()  

    rd = 287.04
    cp = 1004.
    rho_air = 100.*pres/(rd*tv)
    avg_z = np.mean(z,axis=(0,1))

    # Find closest value to 2500 m. AGL
    nearest_val,nearest_id = find_nearest(avg_z,2500.)

    qr_lim = qr[:,:,nearest_id]
    nr_lim = nr[:,:,nearest_id]
    rho_air_lim = rho_air[:,:,nearest_id]

    rho_air_arr.append(rho_air_lim)
    qr_arr.append(qr_lim)
    nr_arr.append(nr_lim)
    time_arr.append(time)


rho_air_arr= np.array(rho_air_arr)
qr_arr = np.array(qr_arr)
nr_arr = np.array(nr_arr)
time_arr = np.array(time_arr)

rain_vars_dict = {'rho_air':rho_air_arr,\
                  'qr':qr_arr,\
                  'nr':nr_arr,\
                  'time_arr':time_arr,\
                  'lat':lat,\
                  'lon':lon}

print('Writing to dictionary...')
savdir = '/glade/scratch/mckenna/AMIE/amie_post/rain_vars/'
f = open(savdir+'baseline_1km_rain_vars_at_2.5-km.p','wb')
pickle.dump(rain_vars_dict,f)
f.close()
print('Completed writing to dictionary. Done.')


