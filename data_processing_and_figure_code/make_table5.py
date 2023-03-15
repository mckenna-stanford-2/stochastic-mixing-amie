#==========================================================================
# Title: make_table5.py
# Author: McKenna W. Stanford
# Utility: Plots strings used to make Table 5, domain- and time-mean
# convective and stratiform rain rates.
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
import pandas as pd
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


path = '/glade/scratch/mckenna/AMIE/amie_post/'

#=========================================================
# 3-km simulations
#=========================================================
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
    
conv_rain_rate_mean_dict = {}
strat_rain_rate_mean_dict = {}
    
print('Beginning 3-km simulations...')
for sim in sim_names:
    print('Simulation:',sim)
    #conv_rain_rate_mean_dict[sim] = []
    #strat_rain_rate_mean_dict[sim] = []
       
    #=========================================================
    # convective and stratiform IDs
    #=========================================================
    pkl_file = open(path+'conv_strat_id/{}_conv_strat_id.p'.format(sim),'rb')
    conv_strat_dict = pickle.load(pkl_file)
    pkl_file.close()

    con3_mask = conv_strat_dict['con3_mask']
    dbz = conv_strat_dict['dbz']

    #=======================================
    # rain rates
    #=======================================

    pkl_file = open(path+'inst_rain_rates/{}_inst_rain_rates.p'.format(sim),'rb')
    tmpfile = pickle.load(pkl_file)
    pkl_file.close()   

    rain_rate = tmpfile['rain_rate']
    time = tmpfile['time']
    nt = len(time)
    conv_rain_rate_mean = []
    strat_rain_rate_mean = []


    for tt in range(nt):
        tmp_rain_rate = rain_rate[tt,:,:]
        tmp_dbz = dbz[tt,:,:]

        tmp_conv_strat_id = con3_mask[tt,:,:]
        conv_id = np.where(tmp_conv_strat_id == 1)
        strat_id = np.where(tmp_conv_strat_id == 2)
        
        if False:
            # Diagnostics
            print(np.shape(conv_id))
            print(np.shape(strat_id))
            print(np.max(tmp_rain_rate))
            print(np.min(tmp_rain_rate))
            print(np.max(tmp_rain_rate[conv_id]))
            print(np.min(tmp_rain_rate[conv_id]))
            print(np.max(tmp_rain_rate[strat_id]))
            print(np.min(tmp_rain_rate[strat_id]))
            print(np.max(tmp_dbz[conv_id]))
            print(np.min(tmp_dbz[conv_id]))
            print(np.max(tmp_dbz[strat_id]))
            print(np.min(tmp_dbz[strat_id]))


        tmp_conv_rain_rate_mean = np.mean(tmp_rain_rate[conv_id])
        tmp_strat_rain_rate_mean = np.mean(tmp_rain_rate[strat_id])

        conv_rain_rate_mean.append(tmp_conv_rain_rate_mean)
        strat_rain_rate_mean.append(tmp_strat_rain_rate_mean)


    conv_rain_rate_mean = np.array(conv_rain_rate_mean)
    strat_rain_rate_mean = np.array(strat_rain_rate_mean)    
    
    conv_rain_rate_mean_dict[sim] = conv_rain_rate_mean
    strat_rain_rate_mean_dict[sim] = strat_rain_rate_mean




# Calculate time means

keys = list(conv_rain_rate_mean_dict.keys())

mean_conv_rain_rate_dict = {}
mean_strat_rain_rate_dict = {}

for key in keys:
    mean_conv_rain_rate_dict[key] = np.mean(conv_rain_rate_mean_dict[key])
    mean_strat_rain_rate_dict[key] = np.mean(strat_rain_rate_mean_dict[key])

    
    
# Calculate relative differences from baseline
rel_diff_conv_rain_rate_dict = {}
rel_diff_strat_rain_rate_dict = {}

for key in keys:
    if key == 'baseline':
        baseline_conv_rain_rate = mean_conv_rain_rate_dict['baseline']
        baseline_strat_rain_rate = mean_strat_rain_rate_dict['baseline']
        continue
        
    rel_diff_conv_rain_rate_dict[key] = {}
    rel_diff_strat_rain_rate_dict[key] = {}
    
    tmp_conv_rain_rate = mean_conv_rain_rate_dict[key]
    tmp_strat_rain_rate = mean_strat_rain_rate_dict[key]

    rel_diff_conv_rain_rate_dict[key] = ((tmp_conv_rain_rate - baseline_conv_rain_rate)/baseline_conv_rain_rate )*100.
    rel_diff_strat_rain_rate_dict[key] = ((tmp_strat_rain_rate - baseline_strat_rain_rate)/baseline_strat_rain_rate )*100.


    
# Make Table 5
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
for key in keys:
    
    if key == 'baseline':
        baseline_conv_rr_str = str(np.around(mean_conv_rain_rate_dict['baseline'],2))
        baseline_strat_rr_str = str(np.around(mean_strat_rain_rate_dict['baseline'],2))
        
        conv_rr_str = baseline_conv_rr_str + amp_str + '-' + amp_str
        strat_rr_str = baseline_strat_rr_str + amp_str + '-' + space_str + slash_str
        
        out_str = sim_strs[dumi]+amp_str+conv_rr_str+strat_rr_str
    else:
        conv_rr_str = str(np.around(mean_conv_rain_rate_dict[key],2))
        strat_rr_str = str(np.around(mean_strat_rain_rate_dict[key],2))
        
        rel_diff_conv_rr_str = str(np.around(rel_diff_conv_rain_rate_dict[key],1))
        if rel_diff_conv_rain_rate_dict[key] > 0.:
            sign_str = '+'
        else:
            sign_str = ''
        rel_diff_conv_rr_str = sign_str + rel_diff_conv_rr_str   
            
        rel_diff_strat_rr_str = str(np.around(rel_diff_strat_rain_rate_dict[key],1))
        if rel_diff_strat_rain_rate_dict[key] > 0.:
            sign_str = '+'
        else:
            sign_str = ''
        rel_diff_strat_rr_str = sign_str + rel_diff_strat_rr_str          
        

        
        conv_rr_str = conv_rr_str + amp_str + rel_diff_conv_rr_str + amp_str
        strat_rr_str = strat_rr_str + amp_str + rel_diff_strat_rr_str + space_str + slash_str
        
        
        out_str = sim_strs[dumi]+amp_str+ conv_rr_str + strat_rr_str

    print(out_str)
        
    dumi+=1

