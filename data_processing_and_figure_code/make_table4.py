#==========================================================================
# Title: make_table4.py
# Author: McKenna W. Stanford
# Utility: Plots strings used to make Table 4, domain- and time-mean
# total echo area, convective area, and stratiform area
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

path = '/glade/scratch/mckenna/AMIE/amie_post/conv_strat_id/'
#tmpfiles = sorted(glob.glob(path+'*.p'))
#files = []
#for tmpfile in tmpfiles:
#    if ('1km_CG' in tmpfile) or ('spol' in tmpfile):
#        continue
#    files.append(tmpfile)
#num_sims = len(files)
#============================

conv_area_dict = {}
strat_area_dict = {}
tot_echo_area_dict = {}

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
             '4x']

num_sims = len(sim_names)

for sim_ii in range(num_sims):
    print('Simulation:',sim_names[sim_ii])
    
    # Load file
    pkl_file = open(path+sim_names[sim_ii]+'_conv_strat_id.p','rb')
    tmpfile = pickle.load(pkl_file)
    pkl_file.close()
    con3_mask = tmpfile['con3_mask'].copy()
    nt = len(con3_mask[:,0,0])
    conv_area = []
    strat_area = []
    tot_echo_area = []
    for tt in range(nt):
        tmp_con3_mask = np.ndarray.flatten(con3_mask[tt,:,:])
        conv_id = np.where(tmp_con3_mask == 1.)
        strat_id = np.where(tmp_con3_mask == 2.)
        tot_echo_id = np.where(tmp_con3_mask > 0.)
        tmp_conv_area = np.size(conv_id)*9.
        tmp_strat_area = np.size(strat_id)*9.
        tmp_tot_echo_area = np.size(tot_echo_id)*9.
        
        conv_area.append(tmp_conv_area)
        strat_area.append(tmp_strat_area)
        tot_echo_area.append(tmp_tot_echo_area)

    conv_area = np.array(conv_area)
    strat_area = np.array(strat_area)
    tot_echo_area = np.array(tot_echo_area)
    
    tot_echo_area_dict[sim_names[sim_ii]] = tot_echo_area
    conv_area_dict[sim_names[sim_ii]] = conv_area
    strat_area_dict[sim_names[sim_ii]] = strat_area
    #print('Tot. Echo area [x10^5 km^2',np.mean(tot_echo_area)*1.e-5)
    #print('Conv. area [x10^5 km^2',np.mean(conv_area)*1.e-5)
    #print('Strat. area [x10^5 km^2',np.mean(strat_area)*1.e-5)
    


        
        
# Plot time series
iplot = True
if iplot:
    time = tmpfile['time']
    time2 = np.array([pd.to_datetime(time[dd]) for dd in range(len(time))])
    time = time2
    sims1 = ['baseline','stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
    sims2 = ['baseline','stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
    sims3 = ['baseline','theta_pert1','theta_pert2','theta_pert3','theta_pert4',\
             'theta_pert5','no_mixing','4x']

    #lws1 = [3,1.5,1.5,1.5,1.5,1.5]
    lws1 = [3,1,1,1,1,1]
    #lws2 = [3,1.5,1.5,1.5,1.5,1.5]
    lws2 = [3,1,1,1,1,1]
    #lws3 = [3,1.5,1.5,1.5,1.5,1.5,3,3]
    lws3 = [3,1,1,1,1,1,3,3]

    colors1 = ['black','lightsalmon','salmon','red','darkred','firebrick']
    colors2 = ['black', 'powderblue','deepskyblue','dodgerblue', 'blue','navy']
    colors3 = ['black','orchid','fuchsia','mediumorchid','purple','darkorchid','black','black'] 

    lss1 = ['solid','solid','solid','solid','solid','solid']
    lss2 = ['solid','solid','solid','solid','solid','solid']
    lss3 = ['solid','solid','solid','solid','solid','solid','dotted','dashed']
    
    
    fig = plt.figure(figsize=(18,12))
    ax1 = fig.add_subplot(331)
    ax2 = fig.add_subplot(332)
    ax3 = fig.add_subplot(333)
    ax4 = fig.add_subplot(334)
    ax5 = fig.add_subplot(335)
    ax6 = fig.add_subplot(336)
    ax7 = fig.add_subplot(337)
    ax8 = fig.add_subplot(338)
    ax9 = fig.add_subplot(339)
    Fontsize=14
    axlist = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]
    
    for ax in axlist:
        ax.set_xlabel('Time',fontsize=Fontsize)
        ax.tick_params(labelsize=Fontsize)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_xticks(time[::48])
        ax.grid()
    ax1.set_ylabel('Total Echo Area [x10$^{5}$ km$^{2}$]',fontsize=Fontsize)
    ax4.set_ylabel('Convective Area [x10$^{5}$ km$^{2}$]',fontsize=Fontsize)
    ax7.set_ylabel('Stratiform Area [x10$^{5}$ km$^{2}$]',fontsize=Fontsize)

    # Total Echo Area
    dum = 0
    for sim in sims1:
        ax1.plot(time,tot_echo_area_dict[sim]*1.e-5,c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
        dum+=1
    dum = 0
    for sim in sims2:
        ax2.plot(time,tot_echo_area_dict[sim]*1.e-5,c=colors2[dum],lw=lws2[dum],ls=lss2[dum])
        dum+=1  
    dum = 0
    for sim in sims3:
        ax3.plot(time,tot_echo_area_dict[sim]*1.e-5,c=colors3[dum],lw=lws3[dum],ls=lss3[dum])
        dum+=1
        
    # Convective Area
    dum = 0
    for sim in sims1:
        ax4.plot(time,conv_area_dict[sim]*1.e-5,c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
        dum+=1
    dum = 0
    for sim in sims2:
        ax5.plot(time,conv_area_dict[sim]*1.e-5,c=colors2[dum],lw=lws2[dum],ls=lss2[dum])
        dum+=1  
    dum = 0
    for sim in sims3:
        ax6.plot(time,conv_area_dict[sim]*1.e-5,c=colors3[dum],lw=lws3[dum],ls=lss3[dum])
        dum+=1  
        
    # Stratiform Area
    dum = 0
    for sim in sims1:
        ax7.plot(time,strat_area_dict[sim]*1.e-5,c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
        dum+=1
    dum = 0
    for sim in sims2:
        ax8.plot(time,strat_area_dict[sim]*1.e-5,c=colors2[dum],lw=lws2[dum],ls=lss2[dum])
        dum+=1  
    dum = 0
    for sim in sims3:
        ax9.plot(time,strat_area_dict[sim]*1.e-5,c=colors3[dum],lw=lws3[dum],ls=lss3[dum])
        dum+=1    
        
        
    plt.subplots_adjust(wspace=0.2)
    plt.show()
    plt.close()


# Calculate time means

keys = list(conv_area_dict.keys())

mean_conv_area_dict = {}
mean_strat_area_dict = {}
mean_tot_echo_area_dict = {}

for key in keys:
    mean_conv_area_dict[key] = np.mean(conv_area_dict[key])
    mean_strat_area_dict[key] = np.mean(strat_area_dict[key])
    mean_tot_echo_area_dict[key] = np.mean(tot_echo_area_dict[key])

    
# Calculate relative differences from baseline
rel_diff_conv_dict = {}
rel_diff_strat_dict = {}
rel_diff_tot_echo_dict = {}

for key in keys:
    if key == 'baseline':
        baseline_conv_area = mean_conv_area_dict['baseline']
        baseline_strat_area = mean_strat_area_dict['baseline']
        baseline_tot_echo_area = mean_tot_echo_area_dict['baseline']
        continue
        
    rel_diff_conv_dict[key] = {}
    rel_diff_strat_dict[key] = {}
    rel_diff_tot_echo_dict[key] = {}
    
    tmp_conv_area = mean_conv_area_dict[key]
    tmp_strat_area = mean_strat_area_dict[key]
    tmp_tot_echo_area = mean_tot_echo_area_dict[key]

    rel_diff_conv_dict[key] = ((tmp_conv_area - baseline_conv_area)/baseline_conv_area )*100.
    rel_diff_strat_dict[key] = ((tmp_strat_area - baseline_strat_area)/baseline_strat_area )*100.
    rel_diff_tot_echo_dict[key] = ((tmp_tot_echo_area - baseline_tot_echo_area)/baseline_tot_echo_area )*100.


    
# Make Table 4
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
        baseline_conv_str = str(np.around(mean_conv_area_dict['baseline']*1.e-5,2))
        baseline_strat_str = str(np.around(mean_strat_area_dict['baseline']*1.e-5,2))
        baseline_tot_echo_str = str(np.around(mean_tot_echo_area_dict['baseline']*1.e-5,2))
        
        tot_echo_str = baseline_tot_echo_str + amp_str + '-' + amp_str
        conv_str = baseline_conv_str + amp_str + '-' + amp_str
        strat_str = baseline_strat_str + amp_str + '-' + space_str + slash_str
        
        out_str = sim_strs[dumi]+amp_str+tot_echo_str+conv_str+strat_str
        
    else:
        conv_str = str(np.around(mean_conv_area_dict[key]*1.e-5,2))
        strat_str = str(np.around(mean_strat_area_dict[key]*1.e-5,2))
        tot_echo_str = str(np.around(mean_tot_echo_area_dict[key]*1.e-5,2))  
        
        rel_diff_conv_str = str(np.around(rel_diff_conv_dict[key],1))
        if rel_diff_conv_dict[key] > 0.:
            sign_str = '+'
        else:
            sign_str = ''
        rel_diff_conv_str = sign_str + rel_diff_conv_str   
            
        rel_diff_strat_str = str(np.around(rel_diff_strat_dict[key],1))
        if rel_diff_strat_dict[key] > 0.:
            sign_str = '+'
        else:
            sign_str = ''
        rel_diff_strat_str = sign_str + rel_diff_strat_str          
        
        rel_diff_tot_echo_str = str(np.around(rel_diff_tot_echo_dict[key],1))
        if rel_diff_tot_echo_dict[key] > 0.:
            sign_str = '+'
        else:
            sign_str = ''
        rel_diff_tot_echo_str = sign_str + rel_diff_tot_echo_str  
        
        
        tot_echo_str = tot_echo_str + amp_str + rel_diff_tot_echo_str + amp_str 
        conv_str = conv_str + amp_str + rel_diff_conv_str + amp_str
        strat_str = strat_str + amp_str + rel_diff_strat_str + space_str + slash_str
        
        
        out_str = sim_strs[dumi]+amp_str+tot_echo_str + conv_str + strat_str

    print(out_str)
        
    dumi+=1

