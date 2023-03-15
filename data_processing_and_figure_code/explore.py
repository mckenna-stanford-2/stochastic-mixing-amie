#==========================================================================
# Title: explore.py
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

sim_names = sim_names[0:6]



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
    conv_up_mf_mean_dict = {}
    conv_up_mf_tot_dict = {}
    conv_area_dict = {}

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
        conv_area = []
        conv_up_mf_mean = []
        conv_up_mf_tot = []

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

            tmp_cond = col_int_cond[tt,:,:]
            tmp_evap = col_int_evap[tt,:,:]
            tmp_pre = pre[tt,:,:]
            #tmp_pe = pe[tt,:,:]

            # Calculate time mean
            tmp_cond_mean = np.nanmean(tmp_cond)
            tmp_evap_mean = np.nanmean(tmp_evap)
            tmp_pre_mean = np.nanmean(tmp_pre)
            #tmp_pe_mean = np.nanmean(tmp_pee)


            conv_id = np.where(tmp_conv_strat_id == 1)
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


            conv_rain_rate_mean.append(tmp_conv_rain_rate_mean)
            conv_area.append(tmp_conv_area)
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
        conv_area = np.array(conv_area) 
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
        conv_area_dict[sim] = conv_area
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
                'conv_area_dict':conv_area_dict,\
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
    f = open(savdir+'stoch_short_time_series.p','wb')
    pickle.dump(out_dict,f)

#==========================================
# To read in 
#==========================================
if iread:
    savdir = '/glade/scratch/mckenna/AMIE/amie_post/stats/'
    pkl_file = open(savdir+'stoch_short_time_series.p','rb')
    out_dict = pickle.load(pkl_file)
    pkl_file.close()   
    
    conv_rr_mean_dict = out_dict['conv_rr_mean_dict']
    conv_area_dict = out_dict['conv_area_dict']
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
# Plot Time Series
#==========================================
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
    
    
    
        
