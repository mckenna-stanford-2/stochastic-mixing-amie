#==========================================================================
# Title: make_table3.py
# Author: McKenna W. Stanford
# Utility: Plots strings used to make Table 3, domain- and time-mean
# vertically integrated process rates
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

path = '/glade/scratch/mckenna/AMIE/amie_post/proc_rates/'


# Load file
pkl_file = open(path+'mean_proc_rates_all_sims.p','rb')
tmpfile = pickle.load(pkl_file)
pkl_file.close()

keys = list(tmpfile['cond_time_mean'].keys())

proc_rate_dict = {}
for key in keys:
    proc_rate_dict[key] = {}
    proc_rate_dict[key]['cond'] = tmpfile['cond_domain_mean'][key]
    proc_rate_dict[key]['evap'] = tmpfile['evap_domain_mean'][key]
    proc_rate_dict[key]['pre'] = tmpfile['pre_domain_mean'][key]
    proc_rate_dict[key]['pe'] = tmpfile['pre_domain_mean'][key]/tmpfile['cond_domain_mean'][key]

    
# Calculate relative differences from baseline
rel_diff_dict = {}
for key in keys:
    if key == 'baseline':
        baseline_cond_val = proc_rate_dict[key]['cond']
        baseline_evap_val = proc_rate_dict[key]['evap']
        baseline_pre_val = proc_rate_dict[key]['pre']
        baseline_pe_val = proc_rate_dict[key]['pe']
        continue
        
    rel_diff_dict[key] = {}
    cond_val = proc_rate_dict[key]['cond']
    evap_val = proc_rate_dict[key]['evap']
    pre_val = proc_rate_dict[key]['pre']
    pe_val = proc_rate_dict[key]['pe']
    
    rel_diff_dict[key]['cond'] = ((cond_val - baseline_cond_val)/baseline_cond_val)*100.
    rel_diff_dict[key]['evap'] = ((evap_val - baseline_evap_val)/baseline_evap_val)*100.
    rel_diff_dict[key]['pre'] = ((pre_val - baseline_pre_val)/baseline_pre_val)*100.
    rel_diff_dict[key]['pe'] = ((pe_val - baseline_pe_val)/baseline_pe_val)*100.

    
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
for key,val in proc_rate_dict.items():
    
    if key == 'baseline':
        cond_str = str(np.around(proc_rate_dict[key]['cond']*1.e4,2))
        cond_str = cond_str + amp_str
        pre_str = str(np.around(proc_rate_dict[key]['pre']*1.e4,2))
        pre_str = pre_str + amp_str
        evap_str = str(np.around(proc_rate_dict[key]['evap']*1.e4,2))
        evap_str = evap_str + amp_str
        pe_str = str(np.around(proc_rate_dict[key]['pe'],2))
        pe_str = pe_str + amp_str + '-' + space_str + slash_str
        out_str = sim_strs[dumi]+amp_str+pre_str+'-'+amp_str+cond_str+'-'+amp_str+evap_str+'-'+amp_str+pe_str
    else:
        cond_str = str(np.around(proc_rate_dict[key]['cond']*1.e4,2))
        evap_str = str(np.around(proc_rate_dict[key]['evap']*1.e4,2))
        pre_str = str(np.around(proc_rate_dict[key]['pre']*1.e4,2))
        pe_str = str(np.around(proc_rate_dict[key]['pe'],2))
        cond_str = cond_str+amp_str
        evap_str = evap_str+amp_str
        pre_str = pre_str+amp_str
        pe_str = pe_str+amp_str
        
        rel_diff_cond_str = str(np.around(rel_diff_dict[key]['cond'],1))
        rel_diff_evap_str = str(np.around(rel_diff_dict[key]['evap'],1))
        rel_diff_pre_str = str(np.around(rel_diff_dict[key]['pre'],1))
        rel_diff_pe_str = str(np.around(rel_diff_dict[key]['pe'],1))
        
        # Pre
        if rel_diff_dict[key]['pre'] > 0.:
            sign_str = '+'
        else:
            sign_str = ''
        rel_diff_pre_str = sign_str + rel_diff_pre_str
        rel_diff_pre_str = rel_diff_pre_str+amp_str
        
        # Cond
        if rel_diff_dict[key]['cond'] > 0.:
            sign_str = '+'
        else:
            sign_str = ''
        rel_diff_cond_str = sign_str + rel_diff_cond_str
        rel_diff_cond_str = rel_diff_cond_str+amp_str        
        
        # Evap
        if rel_diff_dict[key]['evap'] > 0.:
            sign_str = '+'
        else:
            sign_str = ''
        rel_diff_evap_str = sign_str + rel_diff_evap_str
        rel_diff_evap_str = rel_diff_evap_str+amp_str  
        
        # PE
        if rel_diff_dict[key]['pe'] > 0.:
            sign_str = '+'
        else:
            sign_str = ''
        rel_diff_pe_str = sign_str + rel_diff_pe_str
        rel_diff_pe_str = rel_diff_pe_str + space_str + slash_str  
        
        
        out_str = sim_strs[dumi]+amp_str+pre_str+rel_diff_pre_str+cond_str+rel_diff_cond_str+evap_str+rel_diff_evap_str+pe_str+rel_diff_pe_str

    print(out_str)
        
    dumi+=1

