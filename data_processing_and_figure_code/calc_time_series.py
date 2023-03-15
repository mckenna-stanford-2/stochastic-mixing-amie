#==========================================================================
# Title: calc_time_series.py
# Author: McKenna W. Stanford
# Utility: 
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
from scipy import stats as stats
import pandas as pd
from matplotlib import cm
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import PolynomialFeatures
import operator
import pandas as pd
import matplotlib.dates as mdates


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

#------------------------------------------
# Parameters
#------------------------------------------
plt.rc('text',usetex=True)
dfmt = mdates.DateFormatter('%d-%H')


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

#sim_names = sim_names[0:6]

#sim_names = sim_names[0:1]

#---------------------------------------------
#---------------------------------------------
# Option to skip and read in dictionary if
# calculations have already been made
#---------------------------------------------
#---------------------------------------------
icalc = 1
iread = 0

if icalc:
    conv_rr_mean_dict = {}
    strat_rr_mean_dict = {}
    
    up_mf_mean_dict = {}
    conv_up_mf_mean_dict = {}
    conv_up_mf_tot_dict = {}
    
    conv_area_dict = {}
    strat_area_dict = {}
    tot_echo_area_dict = {}

    cond_mean_dict = {}
    evap_mean_dict = {}
    pre_mean_dict = {}
    pe_mean_dict = {}

    conv_cond_mean_dict = {}
    conv_evap_mean_dict = {}
    conv_pre_mean_dict = {}
    conv_pe_mean_dict = {}

    path = '/glade/scratch/mckenna/AMIE/amie_post/'

    for sim in sim_names:
        print('Simulation:',sim)

        #=========================================================
        # convective and stratiform IDs
        #=========================================================
        pkl_file = open(path+'conv_strat_id/{}_conv_strat_id.p'.format(sim),'rb')
        conv_strat_dict = pickle.load(pkl_file)
        pkl_file.close()

        con3_mask = conv_strat_dict['con3_mask']
        con3_mask = con3_mask
        #dbz = conv_strat_dict['dbz']

        #=======================================
        # rain rates
        #=======================================

        pkl_file = open(path+'inst_rain_rates/{}_inst_rain_rates.p'.format(sim),'rb')
        tmpfile = pickle.load(pkl_file)
        pkl_file.close()   

        rain_rate = tmpfile['rain_rate']
        time = tmpfile['time']
        nt = len(time)

        #=======================================
        # Mass Flux
        #=======================================
        pkl_file = open(path+'mass_flux/{}_int_mass_flux_dict.p'.format(sim),'rb')
        int_mf_dict = pickle.load(pkl_file)
        pkl_file.close()
        conv_up_mf = int_mf_dict['conv_up_mf']
        up_mf = int_mf_dict['pos_up_mf']

        #=======================================
        # Proc Rates
        #=======================================
        pkl_file = open(path+'proc_rates/{}_proc_rates.p'.format(sim),'rb')
        proc_rates_dict = pickle.load(pkl_file)
        pkl_file.close()

        col_int_cond = proc_rates_dict['col_int_cond']
        col_int_evap = proc_rates_dict['col_int_evap']
        rainnc = proc_rates_dict['rainnc']
        snownc = proc_rates_dict['snownc']
        time = proc_rates_dict['time']

        pre = rainnc + snownc
        # Convert units
        pre = pre*4./3600.*997./1000. # convert to mass flux
        col_int_cond = col_int_cond*4./3600. # convert to per second
        col_int_evap = col_int_evap*4./3600. # convert to per second  

       # Calculate precipitation efficiency
        pe = np.zeros(np.shape(pre))
        dumid = np.where(col_int_cond > 0.)
        pe[dumid] = pre[dumid]/col_int_cond[dumid]

        # Loop through times

        conv_rain_rate_mean = []
        strat_rain_rate_mean = []
        
        conv_area = []
        strat_area = []
        tot_echo_area = []
        
        conv_up_mf_mean = []
        conv_up_mf_tot = []
        up_mf_mean = []

        cond_mean = []
        evap_mean = []
        pre_mean = []
        pe_mean = []

        conv_cond_mean = []
        conv_evap_mean = []
        conv_pre_mean = []
        conv_pe_mean = []

        #col_int_evap_time_mean = np.nanmean(col_int_evap,axis=(1,2))
        #col_int_cond_time_mean = np.nanmean(col_int_cond,axis=(1,2))
        #pre_time_mean = np.nanmean(pre,axis=(1,2))
        #plt.plot(time,col_int_evap_time_mean)

        for tt in range(nt):
            tmp_rain_rate = rain_rate[tt,:,:]
            tmp_conv_strat_id = con3_mask[tt,:,:]
            tmp_conv_up_mf = conv_up_mf[:,:,tt]
            tmp_up_mf = up_mf[:,:,tt]

            tmp_cond = col_int_cond[tt,:,:]
            tmp_evap = col_int_evap[tt,:,:]
            tmp_pre = pre[tt,:,:]
            #tmp_pe = pe[tt,:,:]

            # Calculate time mean
            tmp_cond_mean = np.nanmean(tmp_cond)
            tmp_evap_mean = np.nanmean(tmp_evap)
            tmp_pre_mean = np.nanmean(tmp_pre)
            #tmp_pe_mean = np.nanmean(tmp_pee)
            tmp_up_mf_mean = np.nanmean(tmp_up_mf)

            conv_id = np.where(tmp_conv_strat_id == 1)
            strat_id = np.where(tmp_conv_strat_id == 2)
            tot_echo_id = np.where(tmp_conv_strat_id > 0.)
            
            # Convective IDs
            if np.size(conv_id) > 0.:
                tmp_conv_rain_rate_mean = np.mean(tmp_rain_rate[conv_id])
                
                tmp_conv_area = np.size(conv_id[0])*9.
                
                tmp_conv_up_mf_mean = np.mean(tmp_conv_up_mf[tmp_conv_up_mf > 0.])
                tmp_conv_up_mf_tot = np.sum(tmp_conv_up_mf)

                tmp_conv_cond_mean = np.mean(tmp_cond[conv_id])
                tmp_conv_evap_mean = np.mean(tmp_evap[conv_id])
                tmp_conv_pre_mean = np.mean(tmp_pre[conv_id])
                #tmp_conv_pe_mean = np.mean(tmp_pe[conv_id])

            else:
                tmp_conv_rain_rate_mean = np.nan
                tmp_conv_area = 0.
                tmp_conv_up_mf_mean = np.nan
                tmp_conv_up_mf_tot = 0.

                tmp_conv_cond_mean = np.nan
                tmp_conv_evap_mean = np.nan
                tmp_conv_pre_mean = np.nan
                #tmp_conv_pe_mean = np.nan

            # Stratiform IDs
            if np.size(strat_id) > 0.:
                tmp_strat_rain_rate_mean = np.mean(tmp_rain_rate[strat_id])                
                tmp_strat_area = np.size(strat_id[0])*9.

            else:
                tmp_strat_rain_rate_mean = np.nan
                tmp_strat_area = 0.
                
            # Total Echo IDs
            if np.size(tot_echo_id) > 0.:
                tmp_tot_echo_area = np.size(tot_echo_id[0])*9.
            else:
                tmp_tot_echo_area = 0.
                
                
            conv_rain_rate_mean.append(tmp_conv_rain_rate_mean)
            strat_rain_rate_mean.append(tmp_strat_rain_rate_mean)
            
            conv_area.append(tmp_conv_area)
            strat_area.append(tmp_strat_area)
            tot_echo_area.append(tmp_tot_echo_area)
            
            up_mf_mean.append(tmp_up_mf_mean)
            
            conv_up_mf_mean.append(tmp_conv_up_mf_mean)
            conv_up_mf_tot.append(tmp_conv_up_mf_tot)

            cond_mean.append(tmp_cond_mean)
            evap_mean.append(tmp_evap_mean)
            pre_mean.append(tmp_pre_mean)
            #pe_mean.append(tmp_pe_mean)

            conv_cond_mean.append(tmp_conv_cond_mean)
            conv_evap_mean.append(tmp_conv_evap_mean)
            conv_pre_mean.append(tmp_conv_pre_mean)
            #conv_pe_mean.append(tmp_conv_pe_mean)

        conv_rain_rate_mean = np.array(conv_rain_rate_mean) 
        strat_rain_rate_mean = np.array(strat_rain_rate_mean) 
        
        conv_area = np.array(conv_area) 
        strat_area = np.array(strat_area) 
        tot_echo_area = np.array(tot_echo_area) 
        
        up_mf_mean = np.array(up_mf_mean)
        
        conv_up_mf_mean = np.array(conv_up_mf_mean)
        conv_up_mf_tot = np.array(conv_up_mf_tot)

        cond_mean = np.array(cond_mean)
        evap_mean = np.array(evap_mean)
        pre_mean = np.array(pre_mean)
        pe_mean = pre_mean/cond_mean
        #pe_mean = np.array(pe_mean)

        conv_cond_mean = np.array(conv_cond_mean)
        conv_evap_mean = np.array(conv_evap_mean)
        conv_pre_mean = np.array(conv_pre_mean)
        conv_pe_mean = conv_pre_mean/conv_cond_mean
        #conv_pe_mean = np.array(conv_pe_mean)
        #plt.plot(conv_pe_mean)

        # save to dictionary
        conv_rr_mean_dict[sim] = conv_rain_rate_mean
        strat_rr_mean_dict[sim] = strat_rain_rate_mean
        
        conv_area_dict[sim] = conv_area
        strat_area_dict[sim] = strat_area
        tot_echo_area_dict[sim] = tot_echo_area
        
        up_mf_mean_dict[sim] = up_mf_mean
        
        conv_up_mf_mean_dict[sim] = conv_up_mf_mean
        conv_up_mf_tot_dict[sim] = conv_up_mf_tot

        cond_mean_dict[sim] = cond_mean
        evap_mean_dict[sim] = evap_mean
        pre_mean_dict[sim] = pre_mean
        pe_mean_dict[sim] = pe_mean

        conv_cond_mean_dict[sim] = conv_cond_mean
        conv_evap_mean_dict[sim] = conv_evap_mean
        conv_pre_mean_dict[sim] = conv_pre_mean
        conv_pe_mean_dict[sim] = conv_pe_mean

    out_dict = {}

    out_dict = {'conv_rr_mean_dict':conv_rr_mean_dict,\
                'strat_rr_mean_dict':strat_rr_mean_dict,\
                'conv_area_dict':conv_area_dict,\
                'strat_area_dict':strat_area_dict,\
                'tot_echo_area_dict':tot_echo_area_dict,\
                'up_mf_mean_dict':up_mf_mean_dict,\
                'conv_up_mf_mean_dict':conv_up_mf_mean_dict,\
                'conv_up_mf_tot_dict':conv_up_mf_tot_dict,\
                'cond_mean_dict':cond_mean_dict,\
                'evap_mean_dict':evap_mean_dict,\
                'pre_mean_dict':pre_mean_dict,\
                'pe_mean_dict':pe_mean_dict,\
                'conv_cond_mean_dict':conv_cond_mean_dict,\
                'conv_evap_mean_dict':conv_evap_mean_dict,\
                'conv_pre_mean_dict':conv_pre_mean_dict,\
                'conv_pe_mean_dict':conv_pe_mean_dict,\

                }

    savdir = '/glade/scratch/mckenna/AMIE/amie_post/stats/'
    f = open(savdir+'time_series_all.p','wb')
    pickle.dump(out_dict,f)

    
    
#==========================================
# To read in 
#==========================================
if iread:
    savdir = '/glade/scratch/mckenna/AMIE/amie_post/stats/'
    pkl_file = open(savdir+'time_series_all.p','rb')
    out_dict = pickle.load(pkl_file)
    pkl_file.close()   
    
    conv_rr_mean_dict = out_dict['conv_rr_mean_dict']
    strat_rr_mean_dict = out_dict['strat_rr_mean_dict']
    
    conv_area_dict = out_dict['conv_area_dict']
    strat_area_dict = out_dict['strat_area_dict']
    tot_echo_area_dict = out_dict['tot_echo_area_dict']
    
    up_mf_mean_dict = out_dict['up_mf_mean_dict']
    
    conv_up_mf_mean_dict = out_dict['conv_up_mf_mean_dict']
    conv_up_mf_tot_dict = out_dict['conv_up_mf_tot_dict']
    
    cond_mean_dict = out_dict['cond_mean_dict']
    evap_mean_dict = out_dict['evap_mean_dict']
    pre_mean_dict = out_dict['pre_mean_dict']
    pe_mean_dict = out_dict['pe_mean_dict']
    conv_cond_mean_dict = out_dict['conv_cond_mean_dict']
    conv_evap_mean_dict = out_dict['conv_evap_mean_dict']
    conv_pre_mean_dict = out_dict['conv_pre_mean_dict']
    conv_pe_mean_dict = out_dict['conv_pe_mean_dict']

    

#==========================================
# Plot process rate scatterplots - v2
#==========================================





fig,_axes = plt.subplots(nrows=2,ncols=4,figsize=(12,6))
axes_flat = np.ndarray.flatten(_axes)

Fontsize=14

for ax in axes_flat:
    ax.tick_params(labelsize=Fontsize)
    ax.xaxis.set_visible(False) # Hide only x axis
    ax.set_xlim(-1,6)


fac = 0.9
_axes[0,0].set_ylabel('PRE [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[0,1].set_ylabel('$\\langle$COND$\\rangle$ [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[0,2].set_ylabel('$\\langle$EVAP$\\rangle$ [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[0,3].set_ylabel('PE',fontsize=Fontsize*fac)

fac = 0.9
_axes[1,0].set_ylabel('Conv. PRE [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[1,1].set_ylabel('Conv. $\\langle$COND$\\rangle$ [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[1,2].set_ylabel('Conv. $\\langle$EVAP$\\rangle$ [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[1,3].set_ylabel('Conv. PE',fontsize=Fontsize*fac)


dumx = 0.025
dumy = 0.9
facx = 1.1
_axes[0,0].text(dumx,dumy,'x10$^{-4}$',transform=_axes[0,0].transAxes,fontsize=Fontsize*facx)
_axes[0,1].text(dumx,dumy,'x10$^{-4}$',transform=_axes[0,1].transAxes,fontsize=Fontsize*facx)
_axes[0,2].text(dumx,dumy,'x10$^{-4}$',transform=_axes[0,2].transAxes,fontsize=Fontsize*facx)
_axes[1,0].text(dumx,dumy,'x10$^{-4}$',transform=_axes[1,0].transAxes,fontsize=Fontsize*facx)
_axes[1,1].text(dumx,dumy,'x10$^{-4}$',transform=_axes[1,1].transAxes,fontsize=Fontsize*facx)
_axes[1,2].text(dumx,dumy,'x10$^{-4}$',transform=_axes[1,2].transAxes,fontsize=Fontsize*facx)

_axes[0,0].set_ylim(0.86,1.15)
_axes[0,1].set_ylim(1.625,2)
_axes[0,2].set_ylim(0.75,1.07)
_axes[0,3].set_ylim(0.46,0.58)
_axes[1,0].set_ylim(25,50)
_axes[1,1].set_ylim(45,75)
_axes[1,2].set_ylim(12.7,14.25)
_axes[1,3].set_ylim(0.55,0.675)

x = np.arange(0,6,1)

#------------------------
# PRE
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(pre_mean_dict[sim])
_axes[0,0].plot(x[0],tmp_mean_dict['baseline']*1.e4,\
                    marker='*',markersize=10,c='k')
_axes[0,0].plot(x[1],tmp_mean_dict['4x']*1.e4,\
                    marker='X',markersize=10,c='k')
_axes[0,0].plot(x[2],tmp_mean_dict['no_mixing']*1.e4,\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
_axes[0,0].plot([x[3],x[3]],[stoch_short_min*1.e4,stoch_short_max*1.e4],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
_axes[0,0].plot([x[4],x[4]],[stoch_long_min*1.e4,stoch_long_max*1.e4],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
_axes[0,0].plot([x[5],x[5]],[theta_pert_min*1.e4,theta_pert_max*1.e4],c='purple',lw=7)
#------------------------
# COND
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(cond_mean_dict[sim])
_axes[0,1].plot(x[0],tmp_mean_dict['baseline']*1.e4,\
                    marker='*',markersize=10,c='k')
_axes[0,1].plot(x[1],tmp_mean_dict['4x']*1.e4,\
                    marker='X',markersize=10,c='k')
_axes[0,1].plot(x[2],tmp_mean_dict['no_mixing']*1.e4,\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
_axes[0,1].plot([x[3],x[3]],[stoch_short_min*1.e4,stoch_short_max*1.e4],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
_axes[0,1].plot([x[4],x[4]],[stoch_long_min*1.e4,stoch_long_max*1.e4],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
_axes[0,1].plot([x[5],x[5]],[theta_pert_min*1.e4,theta_pert_max*1.e4],c='purple',lw=7)
#------------------------
# EVAP
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(evap_mean_dict[sim])
_axes[0,2].plot(x[0],tmp_mean_dict['baseline']*1.e4,\
                    marker='*',markersize=10,c='k')
_axes[0,2].plot(x[1],tmp_mean_dict['4x']*1.e4,\
                    marker='X',markersize=10,c='k')
_axes[0,2].plot(x[2],tmp_mean_dict['no_mixing']*1.e4,\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
_axes[0,2].plot([x[3],x[3]],[stoch_short_min*1.e4,stoch_short_max*1.e4],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
_axes[0,2].plot([x[4],x[4]],[stoch_long_min*1.e4,stoch_long_max*1.e4],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
_axes[0,2].plot([x[5],x[5]],[theta_pert_min*1.e4,theta_pert_max*1.e4],c='purple',lw=7)

#------------------------
# PE
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(pe_mean_dict[sim])
_axes[0,3].plot(x[0],tmp_mean_dict['baseline'],\
                    marker='*',markersize=10,c='k')
_axes[0,3].plot(x[1],tmp_mean_dict['4x'],\
                    marker='X',markersize=10,c='k')
_axes[0,3].plot(x[2],tmp_mean_dict['no_mixing'],\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
_axes[0,3].plot([x[3],x[3]],[stoch_short_min,stoch_short_max],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
_axes[0,3].plot([x[4],x[4]],[stoch_long_min,stoch_long_max],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
_axes[0,3].plot([x[5],x[5]],[theta_pert_min,theta_pert_max],c='purple',lw=7)



#------------------------
# Conv. PRE
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(conv_pre_mean_dict[sim])
_axes[1,0].plot(x[0],tmp_mean_dict['baseline']*1.e4,\
                    marker='*',markersize=10,c='k')
_axes[1,0].plot(x[1],tmp_mean_dict['4x']*1.e4,\
                    marker='X',markersize=10,c='k')
_axes[1,0].plot(x[2],tmp_mean_dict['no_mixing']*1.e4,\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
_axes[1,0].plot([x[3],x[3]],[stoch_short_min*1.e4,stoch_short_max*1.e4],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
_axes[1,0].plot([x[4],x[4]],[stoch_long_min*1.e4,stoch_long_max*1.e4],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
_axes[1,0].plot([x[5],x[5]],[theta_pert_min*1.e4,theta_pert_max*1.e4],c='purple',lw=7)
#------------------------
# Conv. COND
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(conv_cond_mean_dict[sim])
_axes[1,1].plot(x[0],tmp_mean_dict['baseline']*1.e4,\
                    marker='*',markersize=10,c='k')
_axes[1,1].plot(x[1],tmp_mean_dict['4x']*1.e4,\
                    marker='X',markersize=10,c='k')
_axes[1,1].plot(x[2],tmp_mean_dict['no_mixing']*1.e4,\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
_axes[1,1].plot([x[3],x[3]],[stoch_short_min*1.e4,stoch_short_max*1.e4],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
_axes[1,1].plot([x[4],x[4]],[stoch_long_min*1.e4,stoch_long_max*1.e4],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
_axes[1,1].plot([x[5],x[5]],[theta_pert_min*1.e4,theta_pert_max*1.e4],c='purple',lw=7)
#------------------------
# Conv. EVAP
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(conv_evap_mean_dict[sim])
_axes[1,2].plot(x[0],tmp_mean_dict['baseline']*1.e4,\
                    marker='*',markersize=10,c='k')
_axes[1,2].plot(x[1],tmp_mean_dict['4x']*1.e4,\
                    marker='X',markersize=10,c='k')
_axes[1,2].plot(x[2],tmp_mean_dict['no_mixing']*1.e4,\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
_axes[1,2].plot([x[3],x[3]],[stoch_short_min*1.e4,stoch_short_max*1.e4],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
_axes[1,2].plot([x[4],x[4]],[stoch_long_min*1.e4,stoch_long_max*1.e4],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
_axes[1,2].plot([x[5],x[5]],[theta_pert_min*1.e4,theta_pert_max*1.e4],c='purple',lw=7)

#------------------------
# Conv. PE
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(conv_pe_mean_dict[sim])
_axes[1,3].plot(x[0],tmp_mean_dict['baseline'],\
                    marker='*',markersize=10,c='k')
_axes[1,3].plot(x[1],tmp_mean_dict['4x'],\
                    marker='X',markersize=10,c='k')
_axes[1,3].plot(x[2],tmp_mean_dict['no_mixing'],\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
_axes[1,3].plot([x[3],x[3]],[stoch_short_min,stoch_short_max],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
_axes[1,3].plot([x[4],x[4]],[stoch_long_min,stoch_long_max],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
_axes[1,3].plot([x[5],x[5]],[theta_pert_min,theta_pert_max],c='purple',lw=7)

    
    
custom_lines = [Line2D([0],[0],color='k', marker='*', linestyle='None',
                          markersize=dum2, label='BASELINE'),\
                Line2D([0],[0],color='black', marker='X', linestyle='None',
                          markersize=dum2, label='4X'),\
                Line2D([0],[0],color='black', marker='d', linestyle='None',
                          markersize=dum2, label='NO\_MIXING'),\
                Line2D([0],[0],color='red', lw=7,\
                          label='STOCH\_SHORT'),\
                Line2D([0],[0],color='blue', lw=7,\
                          label='STOCH\_LONG'),\
                Line2D([0],[0],color='purple', lw=7,\
                          label='$\\theta$-pert'),\
               ]
_axes[1,1].legend(handles = custom_lines,fontsize=Fontsize,\
                loc='lower center',ncol=3,bbox_to_anchor=(1.25,-0.5))    
    
    
plt.subplots_adjust(wspace=0.45)
plt.show()
plt.close()    



#==========================================
# Plot area and rain rate scatterplots - v2
#==========================================




fig,_axes = plt.subplots(nrows=2,ncols=3,figsize=(12,6))
axes_flat = np.ndarray.flatten(_axes)
fig.delaxes(axes_flat[-1])

Fontsize=16

for ax in axes_flat:
    ax.tick_params(labelsize=Fontsize)
    ax.xaxis.set_visible(False) # Hide only x axis
    ax.set_xlim(-1,6)


fac = 1.
axes_flat[0].set_ylabel('Tot. Echo Area [km$^{2}$]',fontsize=Fontsize*fac)
axes_flat[1].set_ylabel('Conv. Area [km$^{2}$]',fontsize=Fontsize*fac)
axes_flat[2].set_ylabel('Strat. Area [km$^{2}$]',fontsize=Fontsize*fac)
axes_flat[3].set_ylabel('Conv. Rain Rate [mm hr$^{-1}$]',fontsize=Fontsize*fac)
axes_flat[4].set_ylabel('Strat. Rain Rate [mm hr$^{-1}$]',fontsize=Fontsize*fac)



dumx = 0.025
dumy = 0.9
facx = 1.1
axes_flat[0].text(dumx,dumy,'x10$^{5}$',transform=axes_flat[0].transAxes,fontsize=Fontsize*facx)
axes_flat[1].text(dumx,dumy,'x10$^{5}$',transform=axes_flat[1].transAxes,fontsize=Fontsize*facx)
axes_flat[2].text(dumx,dumy,'x10$^{5}$',transform=axes_flat[2].transAxes,fontsize=Fontsize*facx)

x = np.arange(0,6,1)

#------------------------
# Tot Echo Area
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(tot_echo_area_dict[sim])
axes_flat[0].plot(x[0],tmp_mean_dict['baseline']*1.e-5,\
                    marker='*',markersize=10,c='k')
axes_flat[0].plot(x[1],tmp_mean_dict['4x']*1.e-5,\
                    marker='X',markersize=10,c='k')
axes_flat[0].plot(x[2],tmp_mean_dict['no_mixing']*1.e-5,\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
axes_flat[0].plot([x[3],x[3]],[stoch_short_min*1.e-5,stoch_short_max*1.e-5],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
axes_flat[0].plot([x[4],x[4]],[stoch_long_min*1.e-5,stoch_long_max*1.e-5],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
axes_flat[0].plot([x[5],x[5]],[theta_pert_min*1.e-5,theta_pert_max*1.e-5],c='purple',lw=7)

#------------------------
# Conv. Area
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(conv_area_dict[sim])
axes_flat[1].plot(x[0],tmp_mean_dict['baseline']*1.e-5,\
                    marker='*',markersize=10,c='k')
axes_flat[1].plot(x[1],tmp_mean_dict['4x']*1.e-5,\
                    marker='X',markersize=10,c='k')
axes_flat[1].plot(x[2],tmp_mean_dict['no_mixing']*1.e-5,\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
axes_flat[1].plot([x[3],x[3]],[stoch_short_min*1.e-5,stoch_short_max*1.e-5],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
axes_flat[1].plot([x[4],x[4]],[stoch_long_min*1.e-5,stoch_long_max*1.e-5],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
axes_flat[1].plot([x[5],x[5]],[theta_pert_min*1.e-5,theta_pert_max*1.e-5],c='purple',lw=7)

    
    
#------------------------
# Strat. Area
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(strat_area_dict[sim])
axes_flat[2].plot(x[0],tmp_mean_dict['baseline']*1.e-5,\
                    marker='*',markersize=10,c='k')
axes_flat[2].plot(x[1],tmp_mean_dict['4x']*1.e-5,\
                    marker='X',markersize=10,c='k')
axes_flat[2].plot(x[2],tmp_mean_dict['no_mixing']*1.e-5,\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
axes_flat[2].plot([x[3],x[3]],[stoch_short_min*1.e-5,stoch_short_max*1.e-5],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
axes_flat[2].plot([x[4],x[4]],[stoch_long_min*1.e-5,stoch_long_max*1.e-5],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
axes_flat[2].plot([x[5],x[5]],[theta_pert_min*1.e-5,theta_pert_max*1.e-5],c='purple',lw=7)

    
#------------------------
# Conv. RR
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(conv_rr_mean_dict[sim])
axes_flat[3].plot(x[0],tmp_mean_dict['baseline'],\
                    marker='*',markersize=10,c='k')
axes_flat[3].plot(x[1],tmp_mean_dict['4x'],\
                    marker='X',markersize=10,c='k')
axes_flat[3].plot(x[2],tmp_mean_dict['no_mixing'],\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
axes_flat[3].plot([x[3],x[3]],[stoch_short_min,stoch_short_max],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
axes_flat[3].plot([x[4],x[4]],[stoch_long_min,stoch_long_max],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
axes_flat[3].plot([x[5],x[5]],[theta_pert_min,theta_pert_max],c='purple',lw=7)

    
#------------------------
# Strat. RR
#------------------------
tmp_mean_dict = {}
for sim in sims:
    tmp_mean_dict[sim] = np.mean(strat_rr_mean_dict[sim])
axes_flat[4].plot(x[0],tmp_mean_dict['baseline'],\
                    marker='*',markersize=10,c='k')
axes_flat[4].plot(x[1],tmp_mean_dict['4x'],\
                    marker='X',markersize=10,c='k')
axes_flat[4].plot(x[2],tmp_mean_dict['no_mixing'],\
                    marker='d',markersize=10,c='k')
# STOCH_SHORT
stoch_short_arr = []
tmp_sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
for sim in tmp_sims:
    stoch_short_arr.append(tmp_mean_dict[sim])
stoch_short_arr = np.array(stoch_short_arr)
stoch_short_min = np.min(stoch_short_arr)
stoch_short_max = np.max(stoch_short_arr)
axes_flat[4].plot([x[3],x[3]],[stoch_short_min,stoch_short_max],c='red',lw=7)
# STOCH_LONG
stoch_long_arr = []
tmp_sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
for sim in tmp_sims:
    stoch_long_arr.append(tmp_mean_dict[sim])
stoch_long_arr = np.array(stoch_long_arr)
stoch_long_min = np.min(stoch_long_arr)
stoch_long_max = np.max(stoch_long_arr)
axes_flat[4].plot([x[4],x[4]],[stoch_long_min,stoch_long_max],c='blue',lw=7)
# theta-pert
theta_pert_arr = []
tmp_sims = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5']
for sim in tmp_sims:
    theta_pert_arr.append(tmp_mean_dict[sim])
theta_pert_arr = np.array(theta_pert_arr)
theta_pert_min = np.min(theta_pert_arr)
theta_pert_max = np.max(theta_pert_arr)
axes_flat[4].plot([x[5],x[5]],[theta_pert_min,theta_pert_max],c='purple',lw=7)
        
    
    
    
custom_lines = [Line2D([0],[0],color='k', marker='*', linestyle='None',
                          markersize=dum2, label='BASELINE'),\
                Line2D([0],[0],color='black', marker='X', linestyle='None',
                          markersize=dum2, label='4X'),\
                Line2D([0],[0],color='black', marker='d', linestyle='None',
                          markersize=dum2, label='NO\_MIXING'),\
                Line2D([0],[0],color='red', lw=7,\
                          label='STOCH\_SHORT'),\
                Line2D([0],[0],color='blue', lw=7,\
                          label='STOCH\_LONG'),\
                Line2D([0],[0],color='purple', lw=7,\
                          label='$\\theta$-pert'),\
               ]
axes_flat[4].legend(handles = custom_lines,fontsize=Fontsize,\
                loc='right',ncol=1,bbox_to_anchor=(2.35,0.5))     
    
    
plt.subplots_adjust(wspace=0.35)
plt.show()
plt.close()    





#==========================================
# Plot process rate scatterplots
#==========================================

sims1 = ['baseline','stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
sims2 = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
sims3 = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4',\
         'theta_pert5','no_mixing','4x']

sims = sims1 + sims2 + sims3
num_sims = len(sims)

colors1 = ['black','lightsalmon','salmon','red','darkred','firebrick']
colors2 = ['powderblue','deepskyblue','dodgerblue', 'blue','navy']
colors3 = ['orchid','fuchsia','mediumorchid','purple','darkorchid','black','black'] 

markers1 = ['*','^','^','^','^','^']
markers2 = ['v','v','v','v','v']
markers3 = ['o','o','o','o','o','d','X']

dum1 = 8
dum2 = 10
sizes1 = [dum2,dum1,dum1,dum1,dum1,dum1]
sizes2 = [dum1,dum1,dum1,dum1,dum1]
sizes3 = [dum1,dum1,dum1,dum1,dum1,dum2,dum2]

colors = colors1 + colors2 + colors3
markers = markers1 + markers2 + markers3
sizes = sizes1 + sizes2 + sizes3






fig,_axes = plt.subplots(nrows=2,ncols=4,figsize=(12,6))
axes_flat = np.ndarray.flatten(_axes)

Fontsize=14

for ax in axes_flat:
    ax.tick_params(labelsize=Fontsize)
    ax.xaxis.set_visible(False) # Hide only x axis
    ax.set_xlim(-1.5,18.5)


fac = 0.9
_axes[0,0].set_ylabel('PRE [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[0,1].set_ylabel('$\\langle$COND$\\rangle$ [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[0,2].set_ylabel('$\\langle$EVAP$\\rangle$ [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[0,3].set_ylabel('PE',fontsize=Fontsize*fac)

fac = 0.9
_axes[1,0].set_ylabel('Conv. PRE [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[1,1].set_ylabel('Conv. $\\langle$COND$\\rangle$ [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[1,2].set_ylabel('Conv. $\\langle$EVAP$\\rangle$ [kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*fac)
_axes[1,3].set_ylabel('Conv. PE',fontsize=Fontsize*fac)


dumx = 0.025
dumy = 0.9
facx = 1.1
_axes[0,0].text(dumx,dumy,'x10$^{-4}$',transform=_axes[0,0].transAxes,fontsize=Fontsize*facx)
_axes[0,1].text(dumx,dumy,'x10$^{-4}$',transform=_axes[0,1].transAxes,fontsize=Fontsize*facx)
_axes[0,2].text(dumx,dumy,'x10$^{-4}$',transform=_axes[0,2].transAxes,fontsize=Fontsize*facx)
_axes[1,0].text(dumx,dumy,'x10$^{-4}$',transform=_axes[1,0].transAxes,fontsize=Fontsize*facx)
_axes[1,1].text(dumx,dumy,'x10$^{-4}$',transform=_axes[1,1].transAxes,fontsize=Fontsize*facx)
_axes[1,2].text(dumx,dumy,'x10$^{-4}$',transform=_axes[1,2].transAxes,fontsize=Fontsize*facx)

# Set yticks
if False:
    _axes[0,0].set_yticks([0.9,1,1.1])
    _axes[0,1].set_yticks([1.7,1.8,1.9])
    _axes[0,2].set_yticks([0.8,0.9])
    _axes[0,3].set_yticks([0.48,0.5,0.52,0.54,0.56])
    _axes[1,0].set_yticks([30,35,40,45])
    _axes[1,1].set_yticks([50,55,60,65,70])
    _axes[1,2].set_yticks([13,14])
    _axes[1,3].set_yticks([0.6,0.65])

x = np.arange(0,num_sims,1)

keys = list(pre_mean_dict.keys())


for sim_ii in range(num_sims):
    
    #-----------------------
    # Domain Mean
    #-----------------------
    # PRE
    tmp_mean = np.mean(pre_mean_dict[sims[sim_ii]])
    _axes[0,0].plot(x[sim_ii],tmp_mean*1.e4,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])
    # COND
    tmp_mean = np.mean(cond_mean_dict[sims[sim_ii]])
    _axes[0,1].plot(x[sim_ii],tmp_mean*1.e4,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])
    # EVAP
    tmp_mean = np.mean(evap_mean_dict[sims[sim_ii]])
    _axes[0,2].plot(x[sim_ii],tmp_mean*1.e4,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])    

    # PE
    tmp_mean = np.mean(pe_mean_dict[sims[sim_ii]])
    _axes[0,3].plot(x[sim_ii],tmp_mean,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])        
    
    #-----------------------
    # Convective Mean
    #-----------------------
    # PRE
    tmp_mean = np.mean(conv_pre_mean_dict[sims[sim_ii]])
    _axes[1,0].plot(x[sim_ii],tmp_mean*1.e4,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])
    # COND
    tmp_mean = np.mean(conv_cond_mean_dict[sims[sim_ii]])
    _axes[1,1].plot(x[sim_ii],tmp_mean*1.e4,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])
    # EVAP
    tmp_mean = np.mean(conv_evap_mean_dict[sims[sim_ii]])
    _axes[1,2].plot(x[sim_ii],tmp_mean*1.e4,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])    

    # PE
    tmp_mean = np.mean(conv_pe_mean_dict[sims[sim_ii]])
    _axes[1,3].plot(x[sim_ii],tmp_mean,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])   
    
    
    
    
custom_lines = [Line2D([0],[0],color='k', marker='*', linestyle='None',
                          markersize=dum2, label='BASELINE'),\
                Line2D([0],[0],color='red', marker='^', linestyle='None',
                          markersize=dum1, label='STOCH\_SHORT'),\
                Line2D([0],[0],color='blue', marker='v', linestyle='None',
                          markersize=dum1, label='STOCH\_LONG'),\
                Line2D([0],[0],color='purple', marker='o', linestyle='None',
                          markersize=dum1, label='$\\theta$-pert'),\
                Line2D([0],[0],color='black', marker='X', linestyle='None',
                          markersize=dum2, label='4X'),\
                Line2D([0],[0],color='black', marker='d', linestyle='None',
                          markersize=dum2, label='NO\_MIXING'),\
               ]
_axes[1,1].legend(handles = custom_lines,fontsize=Fontsize,\
                loc='lower center',ncol=3,bbox_to_anchor=(1,-0.5))    
    
    
plt.subplots_adjust(wspace=0.45)
plt.show()
plt.close()    

    
#==========================================
# Plot areas and rain rates
#==========================================

sims1 = ['baseline','stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
sims2 = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
sims3 = ['theta_pert1','theta_pert2','theta_pert3','theta_pert4',\
         'theta_pert5','no_mixing','4x']

sims = sims1 + sims2 + sims3
num_sims = len(sims)

colors1 = ['black','lightsalmon','salmon','red','darkred','firebrick']
colors2 = ['powderblue','deepskyblue','dodgerblue', 'blue','navy']
colors3 = ['orchid','fuchsia','mediumorchid','purple','darkorchid','black','black'] 

markers1 = ['*','^','^','^','^','^']
markers2 = ['v','v','v','v','v']
markers3 = ['o','o','o','o','o','d','X']

dum1 = 8
dum2 = 10
sizes1 = [dum2,dum1,dum1,dum1,dum1,dum1]
sizes2 = [dum1,dum1,dum1,dum1,dum1]
sizes3 = [dum1,dum1,dum1,dum1,dum1,dum2,dum2]

colors = colors1 + colors2 + colors3
markers = markers1 + markers2 + markers3
sizes = sizes1 + sizes2 + sizes3




fig,_axes = plt.subplots(nrows=2,ncols=3,figsize=(12,6))
axes_flat = np.ndarray.flatten(_axes)
axes_flat

fig.delaxes(axes_flat[-1])

Fontsize=16

for ax in axes_flat:
    ax.tick_params(labelsize=Fontsize)
    ax.xaxis.set_visible(False) # Hide only x axis
    ax.set_xlim(-1.5,18.5)


fac = 1
axes_flat[0].set_ylabel('Tot. Echo Area [km$^{2}$]',fontsize=Fontsize*fac)
axes_flat[1].set_ylabel('Conv. Area [km$^{2}$]',fontsize=Fontsize*fac)
axes_flat[2].set_ylabel('Strat. Area [km$^{2}$]',fontsize=Fontsize*fac)
axes_flat[3].set_ylabel('Conv. Rain Rate [mm hr$^{-1}$]',fontsize=Fontsize*fac)
axes_flat[4].set_ylabel('Strat. Rain Rate [mm hr$^{-1}$]',fontsize=Fontsize*fac)


dumx = 0.025
dumy = 0.875
facx = 1.1
axes_flat[0].text(dumx,dumy,'x10$^{5}$',transform=axes_flat[0].transAxes,fontsize=Fontsize*facx)
axes_flat[1].text(dumx,dumy,'x10$^{5}$',transform=axes_flat[1].transAxes,fontsize=Fontsize*facx)
axes_flat[2].text(dumx,dumy,'x10$^{5}$',transform=axes_flat[2].transAxes,fontsize=Fontsize*facx)



keys = list(pre_mean_dict.keys())


for sim_ii in range(num_sims):

    # Tot Echo Area
    tmp_mean = np.mean(tot_echo_area_dict[sims[sim_ii]])
    axes_flat[0].plot(x[sim_ii],tmp_mean*1.e-5,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])

    # Conv. Area
    tmp_mean = np.mean(conv_area_dict[sims[sim_ii]])
    axes_flat[1].plot(x[sim_ii],tmp_mean*1.e-5,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])
    
    # Strat. Area
    tmp_mean = np.mean(strat_area_dict[sims[sim_ii]])
    axes_flat[2].plot(x[sim_ii],tmp_mean*1.e-5,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])    
    
    # Conv. rr
    tmp_mean = np.mean(conv_rr_mean_dict[sims[sim_ii]])
    axes_flat[3].plot(x[sim_ii],tmp_mean,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii]) 
    # Strat. rr
    tmp_mean = np.mean(strat_rr_mean_dict[sims[sim_ii]])
    axes_flat[4].plot(x[sim_ii],tmp_mean,\
                    marker=markers[sim_ii],markersize=sizes[sim_ii],\
                    c=colors[sim_ii])     
    
    
custom_lines = [Line2D([0],[0],color='k', marker='*', linestyle='None',
                          markersize=dum2, label='BASELINE'),\
                Line2D([0],[0],color='red', marker='^', linestyle='None',
                          markersize=dum1, label='STOCH\_SHORT'),\
                Line2D([0],[0],color='blue', marker='v', linestyle='None',
                          markersize=dum1, label='STOCH\_LONG'),\
                Line2D([0],[0],color='purple', marker='o', linestyle='None',
                          markersize=dum1, label='$\\theta$-pert'),\
                Line2D([0],[0],color='black', marker='X', linestyle='None',
                          markersize=dum2, label='4X'),\
                Line2D([0],[0],color='black', marker='d', linestyle='None',
                          markersize=dum2, label='NO\_MIXING'),\
               ]
axes_flat[4].legend(handles = custom_lines,fontsize=Fontsize,\
                loc='right',ncol=1,bbox_to_anchor=(2.35,0.5))    
    
    
plt.subplots_adjust(wspace=0.35)
plt.show()
plt.close()    

    
    
    
    
    
    

sims1 = ['baseline','stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
sims2 = ['baseline','stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
sims3 = ['baseline','theta_pert1','theta_pert2','theta_pert3','theta_pert4',\
         'theta_pert5','no_mixing','4x']

lws1 = [3,1,1,1,1,1]
lws2 = [3,1,1,1,1,1]
lws3 = [3,1,1,1,1,1,3,3]

colors1 = ['black','lightsalmon','salmon','red','darkred','firebrick']
colors2 = ['black', 'powderblue','deepskyblue','dodgerblue', 'blue','navy']
colors3 = ['black','orchid','fuchsia','mediumorchid','purple','darkorchid','black','black'] 

lss1 = ['solid','solid','solid','solid','solid','solid']
lss2 = ['solid','solid','solid','solid','solid','solid']
lss3 = ['solid','solid','solid','solid','solid','solid','dotted','dashed']
    
    
    
    

#==========================================
# Plot Time Series of vertically integrated
# process rates
#==========================================

fig,_axes = plt.subplots(nrows=4,ncols=3,figsize=(16,14))
axes_flat = np.ndarray.flatten(_axes)

Fontsize=20

for ax in axes_flat:
    ax.tick_params(labelsize=Fontsize)
    ax.xaxis.set_major_formatter(dfmt)
    ax.set_xticks(time[::48])
    ax.grid(ls='dotted',lw=1,c='dimgrey')
    
_axes[0,0].set_ylabel('Mean PRE\n [x10$^{-4}$ kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*1.25)
_axes[1,0].set_ylabel('Mean $\\langle$COND$\\rangle$\n[x10$^{-4}$ kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*1.25)
_axes[2,0].set_ylabel('Mean $\\langle$EVAP$\\rangle$\n[x10$^{-4}$ kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize*1.25)
_axes[3,0].set_ylabel('Mean PE',fontsize=Fontsize*1.5)

_axes[3,0].set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize*1.25)
_axes[3,1].set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize*1.25)
_axes[3,2].set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize*1.25)

# PRE
dum = 0
for sim in sims1:
    _axes[0,0].plot(time,running_mean(pre_mean_dict[sim],4)*1.e4,lw=lws1[dum],ls=lss1[dum],c=colors1[dum])
    dum+=1
    
dum=0
for sim in sims2:
    _axes[0,1].plot(time,running_mean(pre_mean_dict[sim],4)*1.e4,lw=lws2[dum],ls=lss2[dum],c=colors2[dum])
    dum+=1

dum=0
for sim in sims3:
    _axes[0,2].plot(time,running_mean(pre_mean_dict[sim],4)*1.e4,lw=lws3[dum],ls=lss3[dum],c=colors3[dum])
    dum+=1
    
#COND
dum = 0
for sim in sims1:
    _axes[1,0].plot(time,running_mean(cond_mean_dict[sim],4)*1.e4,lw=lws1[dum],ls=lss1[dum],c=colors1[dum])
    dum+=1
    
dum=0
for sim in sims2:
    _axes[1,1].plot(time,running_mean(cond_mean_dict[sim],4)*1.e4,lw=lws2[dum],ls=lss2[dum],c=colors2[dum])
    dum+=1

dum=0
for sim in sims3:
    _axes[1,2].plot(time,running_mean(cond_mean_dict[sim],4)*1.e4,lw=lws3[dum],ls=lss3[dum],c=colors3[dum])
    dum+=1
    
#EVAP
dum = 0
for sim in sims1:
    _axes[2,0].plot(time,running_mean(evap_mean_dict[sim],4)*1.e4,lw=lws1[dum],ls=lss1[dum],c=colors1[dum])
    dum+=1
    
dum=0
for sim in sims2:
    _axes[2,1].plot(time,running_mean(evap_mean_dict[sim],4)*1.e4,lw=lws2[dum],ls=lss2[dum],c=colors2[dum])
    dum+=1

dum=0
for sim in sims3:
    _axes[2,2].plot(time,running_mean(evap_mean_dict[sim],4)*1.e4,lw=lws3[dum],ls=lss3[dum],c=colors3[dum])
    dum+=1    
    
#PE
dum = 0
for sim in sims1:
    _axes[3,0].plot(time,running_mean(pre_mean_dict[sim]/cond_mean_dict[sim],4),\
                    lw=lws1[dum],ls=lss1[dum],c=colors1[dum])
    dum+=1
    
dum=0
for sim in sims2:
    _axes[3,1].plot(time,running_mean(pre_mean_dict[sim]/cond_mean_dict[sim],4),\
                    lw=lws2[dum],ls=lss2[dum],c=colors2[dum])
    dum+=1

dum=0
for sim in sims3:
    _axes[3,2].plot(time,running_mean(pre_mean_dict[sim]/cond_mean_dict[sim],4)\
                    ,lw=lws3[dum],ls=lss3[dum],c=colors3[dum])
    dum+=1    
        
    
dum_max = 3
dum_min= 0.45
_axes[0,0].set_ylim(dum_min,dum_max)
_axes[0,1].set_ylim(dum_min,dum_max)
_axes[0,2].set_ylim(dum_min,dum_max)
_axes[1,0].set_ylim(dum_min,dum_max)
_axes[1,1].set_ylim(dum_min,dum_max)
_axes[1,2].set_ylim(dum_min,dum_max)
_axes[2,0].set_ylim(dum_min,dum_max)
_axes[2,1].set_ylim(dum_min,dum_max)
_axes[2,2].set_ylim(dum_min,dum_max)

dum_min = 0.425
dum_max = 0.62
_axes[3,0].set_ylim(dum_min,dum_max)
_axes[3,1].set_ylim(dum_min,dum_max)
_axes[3,2].set_ylim(dum_min,dum_max)

plt.subplots_adjust(hspace=0.25,wspace=0.25)


plt.show()
plt.close()    
    







#==========================================
# Plot Time Series
#==========================================


    
fig = plt.figure(figsize=(16,12))
ax1 = fig.add_subplot(431)
ax2 = fig.add_subplot(434)
ax3 = fig.add_subplot(437)
ax4 = fig.add_subplot(4,3,10)
ax5 = fig.add_subplot(432)
ax6 = fig.add_subplot(435)
ax7 = fig.add_subplot(438)
ax8 = fig.add_subplot(4,3,11)
ax9 = fig.add_subplot(433)
ax10 = fig.add_subplot(4,3,6)
ax11 = fig.add_subplot(4,3,9)
ax12 = fig.add_subplot(4,3,12)

Fontsize=14
axlist = [ax1,ax2,ax3,ax4,ax5,ax6,\
         ax7,ax8,ax9,ax10,ax11,ax12]

for ax in axlist:
    ax.set_xlabel('Time',fontsize=Fontsize)
    ax.tick_params(labelsize=Fontsize)
    ax.xaxis.set_major_formatter(dfmt)
    ax.set_xticks(time[::48])
    ax.grid()
    
ax1.set_title('Mean Convective Rain Rate',fontsize=Fontsize*1.25)
ax1.set_ylabel('Mean Convective\nRain Rate [mm hr$^{-1}$]',fontsize=Fontsize)

ax2.set_title('Convective Area',fontsize=Fontsize*1.25)
ax2.set_ylabel('Convective Area [km$^{2}$]',fontsize=Fontsize)

ax3.set_title('Mean Convective $\\langle$MF$\\rangle$',fontsize=Fontsize*1.25)
ax3.set_ylabel('Mean Convective\n$\\langle$MF$\\rangle$ [kg s$^{-1}$ m$^{-2}$]',fontsize=Fontsize)

ax4.set_title('Total Convective $\\langle$MF$\\rangle$',fontsize=Fontsize*1.25)
ax4.set_ylabel('Total Convective\n$\\langle$MF$\\rangle$ [kg s$^{-1}$ m$^{-2}$]',fontsize=Fontsize)

ax5.set_title('Mean PRE',fontsize=Fontsize*1.25)
ax5.set_ylabel('Mean PRE [kg s$^{-1}$ m$^{-2}$]',fontsize=Fontsize)

ax6.set_title('Mean $\\langle$COND$\\rangle$',fontsize=Fontsize*1.25)
ax6.set_ylabel('Mean $\\langle$COND$\\rangle$ [kg s$^{-1}$ m$^{-2}$]',fontsize=Fontsize)

ax7.set_title('Mean $\\langle$EVAP$\\rangle$',fontsize=Fontsize*1.25)
ax7.set_ylabel('Mean $\\langle$EVAP$\\rangle$ [kg s$^{-1}$ m$^{-2}$]',fontsize=Fontsize)


ax8.set_title('Mean PE',fontsize=Fontsize*1.25)
ax8.set_ylabel('Mean PE',fontsize=Fontsize)

ax9.set_title('Mean Conv. PRE',fontsize=Fontsize*1.25)
ax9.set_ylabel('Mean Conv.\nPRE [kg s$^{-1}$ m$^{-2}$]',fontsize=Fontsize)

ax10.set_title('Mean Conv. $\\langle$COND$\\rangle$',fontsize=Fontsize*1.25)
ax10.set_ylabel('Mean Conv.\n$\\langle$COND$\\rangle$ [kg s$^{-1}$ m$^{-2}$]',fontsize=Fontsize)

ax11.set_title('Mean Conv. $\\langle$EVAP$\\rangle$',fontsize=Fontsize*1.25)
ax11.set_ylabel('Mean Conv.\n$\\langle$EVAP$\\rangle$ [kg s$^{-1}$ m$^{-2}$]',fontsize=Fontsize)

ax12.set_title('Mean Conv. PE',fontsize=Fontsize*1.25)
ax12.set_ylabel('Mean Conv. PE',fontsize=Fontsize)

dum = 0
for sim in sims1:
    ax1.plot(time,running_mean(conv_rr_mean_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax2.plot(time,running_mean(conv_area_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax3.plot(time,running_mean(conv_up_mf_mean_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax4.plot(time,running_mean(conv_up_mf_tot_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax5.plot(time,running_mean(pre_mean_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax6.plot(time,running_mean(cond_mean_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax7.plot(time,running_mean(evap_mean_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax8.plot(time,running_mean(pe_mean_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax9.plot(time,running_mean(conv_pre_mean_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax10.plot(time,running_mean(conv_cond_mean_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax11.plot(time,running_mean(conv_evap_mean_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax12.plot(time,running_mean(conv_pe_mean_dict[sim],4),c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    #ax4.plot(time,conv_up_mf_mean_dict[sim]*conv_area_dict[sim],c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    dum+=1


plt.subplots_adjust(hspace=0.5,wspace=0.4)
plt.show()
plt.close()    
    
    
    
        
