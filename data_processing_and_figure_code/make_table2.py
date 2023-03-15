#==========================================================================
# Title: make_table2.py
# Author: McKenna W. Stanford
# Utility: Plots strings used to make Table 2, which is domain-accumulated
# volumetric precipitation.
#==========================================================================


#------------------------------------------
# Imports
#------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import glob
import xarray
import matplotlib
import pickle
from netCDF4 import Dataset
import datetime
import wrf
from scipy.ndimage import label
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
#------------------------------------------
# Parameters
#------------------------------------------
dfmt = mdates.DateFormatter('%d-%H')
plt.rc('text',usetex=True)


#------------------------------------------
# Functions
#------------------------------------------
# make function that finds the nearest
# element in array to the given value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx
# Running Mean
def running_mean(x,N):
    x_padded = np.pad(x, (N//2, N-1-N//2), mode='edge')
    x_smooth = np.convolve(x_padded, np.ones((N,))/N, mode='valid')
    return x_smooth

path = '/glade/scratch/mckenna/AMIE/amie_post/acc_gridscale_precip/'


#=========================================================
# Make dictionary to hold all time series
vol_acc_dict = {}
#=========================================================


#=========================================================
# Calculate domain-accumulated volumetric
# precipitation for 3-km simulations.
# Optional in case you've already
# calculated it.
#=========================================================
icalc = 1.
if icalc == 1.:
    print('Processing 3-km simulations...')
    
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

    for sim in sim_names:
        print('Simulation:',sim)
    
        # Load file
        pkl_file = open(path+sim+'_acc_gridscale_precip.p','rb')
        tmpfile = pickle.load(pkl_file)
        pkl_file.close()      

        tmp_rain_rate = tmpfile['rain_rate'][48:,:,:]
        
        # Sum across domain
        domain_sum = np.sum(tmp_rain_rate,axis=(1,2))
        domain_cumsum = np.cumsum(domain_sum)
        # Subtract 1st time since that is what's accumulated during spin-up
        #domain_cumsum = domain_cumsum - domain_cumsum[0]
        # Volumetric precipitatio (mm km^2)
        rain_rate_sum = domain_cumsum[-1]*9.
        vol_acc_dict[sim] = rain_rate_sum

            
        
print('Diagnostics')   
for key,val in vol_acc_dict.items():
    print(key,val*1.e-6)
        
# Calculate relative differences from baseline
rel_diff_dict = {}
for key,val in vol_acc_dict.items():
    if key == 'baseline':
        continue
    rel_diff_dict[key] = ((val - vol_acc_dict['baseline'])/vol_acc_dict['baseline'])*100.
    
for key,val in rel_diff_dict.items():
    print(key,val)
    
# Make Table 1
amp_str = ' & '
slash_str = '\\\\'
space_str = ' '

sim_strs = ['BASELINE',\
           'STOCH1\_SHORT',\
           'STOCH2\_SHORT',\
           'STOCH3\_SHORT',\
           'STOCH4\_SHORT',\
           'STOCH5\_SHORT',\
           'STOCH1\_LONG',\
           'STOCH2\_LONG',\
           'STOCH3\_LONG',\
           'STOCH4\_LONG',\
           'STOCH5\_LONG',\
           '$\\theta$-pert1',\
           '$\\theta$-pert2',\
           '$\\theta$-pert3',\
           '$\\theta$-pert4',\
           '$\\theta$-pert5',\
           'NO\_MIXING',\
           '4X',\
           ]

dumi = 0
for key,val in vol_acc_dict.items():
    
    if key == 'baseline':
        vol_acc_str = str(np.around(vol_acc_dict[key]*1.e-6,2))
        out_str = sim_strs[dumi]+amp_str+vol_acc_str+amp_str+'-'+space_str+slash_str
    else:
        vol_acc_str = str(np.around(vol_acc_dict[key]*1.e-6,2))
        rel_diff_str = str(np.around(rel_diff_dict[key],1))
        if rel_diff_dict[key] > 0.:
            sign_str = '+'
        else:
            sign_str = ''
        out_str = sim_strs[dumi]+amp_str+vol_acc_str+amp_str+sign_str+rel_diff_str+space_str+slash_str
        
    print(out_str)
        
    dumi+=1

