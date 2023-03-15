#==========================================================================
# Title: compute_wrf_acc_gridscale_precip_1km.py
# Author: McKenna W. Stanford
# Utility: Reads in wrfout_d02* (1km sim) files and grabs the total accumulated total
# grid-scale precip (RAINNC) (in mm) and subtracts it from the previous file
# in time, such that the output is accumulated precip per 15-minutes (outfile
# time step). Writes this to a dictionary named 'baseline_1km_acc_gridscale_precip.p'.
#==========================================================================

#------------------------------------------
# Imports
#------------------------------------------
import numpy as np
import glob
import xarray
import pickle
import pandas as pd
import matplotlib.pyplot as plt

#=========================================================
# WRF accumulated total grid-scale precipitation
#=========================================================

path = '/glade/scratch/mckenna/AMIE/baseline_1km/'
files = sorted(glob.glob(path+'wrfout*'))
nt = len(files) # number of files/times

rainnc_arr = []
time_arr = []
for tt in range(nt):
    print('Time step {}/{}:'.format(tt+1,nt))
    ncfile = xarray.open_dataset(files[tt])
    rainnc = ncfile['RAINNC'].values.squeeze().T
    if tt == 0.:
        lat = ncfile['XLAT'].values.squeeze().T
        lon = ncfile['XLONG'].values.squeeze().T  
        lon_1d = lon[:,0]
        lat_1d = lat[0,:]
    time = pd.to_datetime(ncfile['XTIME'].values[0])
    ncfile.close()
    time_arr.append(time)
    rainnc_arr.append(rainnc)

rainnc_arr = np.array(rainnc_arr)
time_arr = np.array(time_arr)

# RAINNC is accumualted in time for each output file
# Need to subtract to get values per time dimension
# So, these will have values of mm/(15 minutes)
rainnc_arr_sub = []
for tt in range(nt-1):
    rainnc_arr_sub.append(rainnc_arr[tt+1,:,:]-rainnc_arr[tt,:,:])
rainnc_arr_sub = np.array(rainnc_arr_sub)

# Limit time to one time step past the start
time_arr = time_arr[1:]

#plt.plot(time_arr,np.mean(rainnc_arr_sub,axis=(1,2)))


# Put in dictionary and write to pickle file
out_dict = {'rain_rate':rainnc_arr_sub,\
            'time':time_arr,\
            'lat':lat,\
            'lon':lon,\
            'lat_1d':lat_1d,\
            'lon_1d':lon_1d}

savdir = '/glade/scratch/mckenna/AMIE/amie_post/acc_gridscale_precip/'
f = open(savdir+'baesline_1km_acc_gridscale_precip.p','wb')
pickle.dump(out_dict,f)
f.close()

    
    
    
    
             



