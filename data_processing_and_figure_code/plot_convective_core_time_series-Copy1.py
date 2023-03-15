#==========================================================================
# Title: plot_convective_cell_stats.py
# Author: McKenna W. Stanford
# Origin Date: 03/20/2020
# Date Modified: 04/28/2020
# Utility: Plots statitsics of observed and simulated convective cells.
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
import pickle

x0 = np.arange(-148.5,151.5,3)
y0 = np.arange(-148.5,151.5,3)
X0,Y0 = np.meshgrid(x0,y0)
radial_dist = np.sqrt(X0**2. + Y0**2.)

savdir = '/glade/scratch/mckenna/AMIE/amie_post/'

#===============================================
# Calculate domain mean multiplicative factor
#===============================================
stochrun = 0
if stochrun == 1:
    median_mult_fac_dict = {}
    mean_mult_fac_dict = {}
    sims = ['stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short',\
           'stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
    
    for sim in sims:
        

        pkl_file = open(savdir+sim+'_convective_objects','rb')
        wrf_objects = pickle.load(pkl_file)
        pkl_file.close() 
    
        tmp = wrf_objects['rand_pert'][:,:,:,0]
        tmpid = np.where(radial_dist > 150.)
        tmp[:,tmpid[0],tmpid[1]] = np.nan
        mult_fac = 2.**(tmp)
        mean_mult_fac = np.nanmean(mult_fac,axis=(1,2))
        median_mult_fac = np.nanmedian(mult_fac,axis=(1,2))
        
        median_mult_fac_dict[sim] = median_mult_fac
        mean_mult_fac_dict[sim] = mean_mult_fac

xtmprun = 0
if xtmprun == 1:
    savdir = '/glade/scratch/mckenna/AMIE/amie_post/'
    pkl_file = open(savdir+'convective_objects_stats_baseline_1km_CG','rb')
    convective_objects_baseline_1km_CG = pickle.load(pkl_file)
    pkl_file.close()
    
tmprun = 0
if tmprun == 1:
    savdir = '/glade/scratch/mckenna/AMIE/amie_post/'
    pkl_file = open(savdir+'convective_objects_stats','rb')
    convective_objects = pickle.load(pkl_file)
    pkl_file.close()
    
    #convective_objects['baseline_1km_CG'] = convective_objects_baseline_1km_CG

    for key,val in convective_objects.items():
        print(key)
        #for key2,val2 in val.items():
        #    print(key2,np.shape(val2))
        #print(' ')


    area = convective_objects['area']
    time = convective_objects['time']
    #if key != 'SPOL':
        #max_w = convective_objects['max_w']
        #avg_w = convective_objects['avg_w']
        #avg_mf = convective_objects['avg_mf']
        #avg_tr = convective_objects['avg_tr']
        #avg_cconevp = convective_objects['avg_cconevp']
        #avg_hdi = convective_objects['avg_hdi']
        #avg_twc = convective_objects['avg_twc']
        #max_twc = convective_objects['max_twc']
        #avg_rand_pert = convective_objects['avg_rand_pert']
    time = convective_objects['time']
    avg_rain_rate = convective_objects['avg_rain_rate']
    max_rain_rate = convective_objects['max_rain_rate']
    #max_echo_top_0 = convective_objects['max_echo_top_0']
    #max_echo_top_10 = convective_objects['max_echo_top_10']
    #max_echo_top_20 = convective_objects['max_echo_top_20']
    #max_echo_top_40 = convective_objects['max_echo_top_40']
    
    area['baseline_1km_CG'] = convective_objects_baseline_1km_CG['area']['baseline_1km_CG']
    #max_w['baseline_1km_CG'] = convective_objects_baseline_1km_CG['max_w']['baseline_1km_CG']
    #avg_w['baseline_1km_CG'] = convective_objects_baseline_1km_CG['avg_w']['baseline_1km_CG']
    #avg_mf['baseline_1km_CG'] = convective_objects_baseline_1km_CG['avg_mf']['baseline_1km_CG']
    #avg_tr['baseline_1km_CG'] = convective_objects_baseline_1km_CG['avg_tr']['baseline_1km_CG']
    #avg_cconevp['baseline_1km_CG'] = convective_objects_baseline_1km_CG['avg_cconevp']['baseline_1km_CG']
    #avg_hdi['baseline_1km_CG'] = convective_objects_baseline_1km_CG['avg_hdi']['baseline_1km_CG']
    #avg_twc['baseline_1km_CG'] = convective_objects_baseline_1km_CG['avg_twc']['baseline_1km_CG']
    #max_twc['baseline_1km_CG'] = convective_objects_baseline_1km_CG['max_twc']['baseline_1km_CG']
    time['baseline_1km_CG'] = convective_objects_baseline_1km_CG['time']['baseline_1km_CG']
    avg_rain_rate['baseline_1km_CG'] = convective_objects_baseline_1km_CG['avg_rain_rate']['baseline_1km_CG']
    max_rain_rate['baseline_1km_CG'] = convective_objects_baseline_1km_CG['max_rain_rate']['baseline_1km_CG']
    #max_echo_top_0['baseline_1km_CG'] = convective_objects_baseline_1km_CG['max_echo_top_0']['baseline_1km_CG']
    #max_echo_top_10['baseline_1km_CG'] = convective_objects_baseline_1km_CG['max_echo_top_10']['baseline_1km_CG']
    #max_echo_top_20['baseline_1km_CG'] = convective_objects_baseline_1km_CG['max_echo_top_20']['baseline_1km_CG']
    #max_echo_top_40['baseline_1km_CG'] = convective_objects_baseline_1km_CG['max_echo_top_40']['baseline_1km_CG']
    #print(aaaaa)
    
    #filter times# let's start with the 7th
    for key,val in time.items():
        print(key)
        val = np.array(val)
        dumval = np.sort(np.unique(val))
        dum7 = dumval[0]
        dum8 = dumval[48]
        dum721 = dumval[37]
        
        #tmpid = np.where(val <= dum721)    
        #tmpid = np.where(val < dum8)
        #tmpid = np.where(val >= dum8)
        tmpid = np.where(val >= dum7)
        
        area[key] = np.array(area[key])[tmpid]
        avg_rain_rate[key] = np.array(avg_rain_rate[key])[tmpid]
        max_rain_rate[key] = np.array(max_rain_rate[key])[tmpid]
        #max_echo_top_0[key] = np.array(max_echo_top_0[key])[tmpid]
        #max_echo_top_10[key] = np.array(max_echo_top_10[key])[tmpid]
        #max_echo_top_20[key] = np.array(max_echo_top_20[key])[tmpid]
        #max_echo_top_40[key] = np.array(max_echo_top_40[key])[tmpid]
        #if key != 'SPOL':
        #    max_w[key] = np.array(max_w[key])[tmpid]
        #    avg_w[key] = np.array(avg_w[key])[tmpid]
        #    avg_mf[key] = np.array(avg_mf[key])[tmpid]
        #    avg_tr[key] = np.array(avg_tr[key])[tmpid]
        #    avg_twc[key] = np.array(avg_twc[key])[tmpid]
        #    max_twc[key] = np.array(max_twc[key])[tmpid]
        #    avg_hdi[key] = np.array(avg_hdi[key])[tmpid]*60. # Convert from K/s to K/min.
        #    avg_cconevp[key] = np.array(avg_cconevp[key])[tmpid]/1005.*60.#Convert from J/K/s to K/min.      
    
    
    # filter sizes
    ifilter2 = 0
    if ifilter2 == 1:
        for key in area.keys():
            print(key)
            tmp_area = np.array(area[key])
            tmp_time = np.array(time[key])
            tmp_avg_rain_rate = np.array(avg_rain_rate[key])
            tmp_max_rain_rate = np.array(max_rain_rate[key])
            tmp_max_echo_top_0 = np.array(max_echo_top_0[key])
            tmp_max_echo_top_20 = np.array(max_echo_top_20[key])
            tmp_max_echo_top_40 = np.array(max_echo_top_40[key])

            tmpid = np.where(tmp_area > 9.)
            tmp_area = tmp_area[tmpid]
            tmp_avg_rain_rate = tmp_avg_rain_rate[tmpid]
            tmp_max_rain_rate = tmp_max_rain_rate[tmpid]
            tmp_max_echo_top_0 = tmp_max_echo_top_0[tmpid]
            tmp_max_echo_top_20 = tmp_max_echo_top_20[tmpid]
            tmp_max_echo_top_40 = tmp_max_echo_top_40[tmpid]
            tmp_time = tmp_time[tmpid]

            area[key] = tmp_area
            time[key] = tmp_time
            avg_rain_rate[key] = tmp_avg_rain_rate
            max_rain_rate[key] = tmp_max_rain_rate
            max_echo_top_0[key] = tmp_max_echo_top_0
            max_echo_top_20[key] = tmp_max_echo_top_20
            max_echo_top_40[key] = tmp_max_echo_top_40
            #print(aaaaa)



#=============================================================
#=============================================================
#=============================================================
#=============================================================
#=============================================================
# Plots
#=============================================================
#=============================================================
#=============================================================
#=============================================================
#=============================================================

colors = ['black',\
         'lightsalmon',\
         'salmon',\
         'red',\
         'darkred',\
         'firebrick',\
         'powderblue',\
         'deepskyblue',\
         'dodgerblue',\
         'blue',\
         'navy',\
         'orchid',\
         'fuchsia',\
         'mediumorchid',\
         'purple',\
         'darkorchid',\
          'green',\
          'lime',\
          'gold',\
          'darkorange',\
         ]  

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
         'baseline_1km_CG',\
            ]
#=========================================
#running mean
#=========================================
def running_mean(x,N):
    x_padded = np.pad(x, (N//2, N-1-N//2), mode='edge')
    x_smooth = np.convolve(x_padded, np.ones((N,))/N, mode='valid')
    return x_smooth
N = 4      

dum_area = {}
dum_time = {}
def_uniq_time = np.unique(time['baseline'])
def_uniq_time_sorted = sorted(def_time)

for sim in sims3:
    tmp_time = convective_objects['time'][sim]
    tmp_time = np.array(tmp_time)
    tmp_area = area[sim]
    tmp_area = np.array(tmp_area)

    #uniq_time = np.unique(tmp_time)
    #uniq_time_sorted = sorted(uniq_time)
    nt = len(def_uniq_time_sorted)
    mean_area_time_series = []
    for tt in range(nt):
        tmpid = np.where(tmp_time == def_uniq_time_sorted[tt])
        if np.size(tmpid) == 0.:
            mean_area_time_series.append(np.nan)
        else:
            tmp_mean_area = np.nanmean(tmp_area[tmpid[0]])
            mean_area_time_series.append(tmp_mean_area)
    mean_area_time_series = np.array(mean_area_time_series)
    dum_area[sim] = mean_area_time_series

    
dum_diff = (dum_area['baseline']-dum_area['no_mixing'])/(dum_area['no_mixing'])*100.
#dum_diff = (dum_area['baseline']-dum_area['no_mixing'])
plt.plot(def_uniq_time_sorted,dum_diff)
    
if True:
    #----------------------------------
    #----------------------------------
    #----------------------------------
    #----------------------------------
    # Time Series of Core Number
    #----------------------------------
    #----------------------------------
    #----------------------------------
    #----------------------------------
    import matplotlib.dates as mdates
    plt.rc('text',usetex=True)
    dfmt = mdates.DateFormatter('%d-%H')

    sims1 = ['baseline',\
             'stoch1_short',\
             'stoch2_short',\
             'stoch3_short',\
             'stoch4_short',\
             'stoch5_short',\
             'SPOL',\
            ]

    sims4 = ['stoch1_short',\
             'stoch2_short',\
             'stoch3_short',\
             'stoch4_short',\
             'stoch5_short']

    sims2 = ['baseline',\
             'stoch1_long',\
             'stoch2_long',\
             'stoch3_long',\
             'stoch4_long',\
             'stoch5_long',\
             'SPOL']

    sims5 = ['stoch1_long',\
             'stoch2_long',\
             'stoch3_long',\
             'stoch4_long',\
             'stoch5_long']         

    sims3 = ['baseline',\
             'theta_pert1',\
             'theta_pert2',\
             'theta_pert3',\
             'theta_pert4',\
             'theta_pert5',\
             'SPOL',\
             'no_mixing',\
             '4x',\
             'baseline_1km_CG',\
            ]

    colors1 = ['black',\
             'lightsalmon',\
             'salmon',\
             'red',\
             'darkred',\
             'firebrick',\
              'grey']

    colors4 = ['lightsalmon',\
             'salmon',\
             'red',\
             'darkred',\
             'firebrick']

    colors2 = ['black',\
             'powderblue',\
             'deepskyblue',\
             'dodgerblue',\
             'blue',\
             'navy',\
             'grey',\
             ]  

    colors5 = ['powderblue',\
             'deepskyblue',\
             'dodgerblue',\
             'blue',\
             'navy'] 

    colors3 = ['black',\
             'orchid',\
             'fuchsia',\
             'mediumorchid',\
             'purple',\
             'darkorchid',\
             'grey',\
             'black',\
             'black',\
             'darkorange',\
             ]  

    lws1 = [4,2,2,2,2,2,4]
    lws2 = [4,2,2,2,2,2,4]
    lws3 = [4,2,2,2,2,2,4,4,4,4]
    lss3 = ['solid','solid','solid','solid','solid','solid','solid','dotted','dashed','-.']

    fig = plt.figure(figsize=(18,14))

    ax1 = fig.add_subplot(331)
    ax2 = fig.add_subplot(332)
    ax3 = fig.add_subplot(333)
    ax4 = fig.add_subplot(334)
    ax5 = fig.add_subplot(335)
    ax6 = fig.add_subplot(336)
    ax7 = fig.add_subplot(337)
    ax8 = fig.add_subplot(338)
    ax9 = fig.add_subplot(339)
    #ax10 = fig.add_subplot(4,3,10)
    #ax11 = fig.add_subplot(4,3,11)

    #axlist = [ax1,ax2,ax3,ax4,ax5,ax6]

    #===============================
    # Number of Cores
    #===============================
    # STOCH_SHORT
    dum = 0
    for sim in sims1:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        core_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            core_time_series.append(np.size(tmpid[0]))
        core_time_series = np.array(core_time_series)
        tmp_core_time_series = running_mean(core_time_series,N)
        ax1.plot(uniq_time_sorted,tmp_core_time_series,lw=lws1[dum],ls='solid',color=colors1[dum])
        dum+=1

    # STOCH_LONG    
    dum = 0
    for sim in sims2:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        core_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            core_time_series.append(np.size(tmpid[0]))
        core_time_series = np.array(core_time_series)
        tmp_core_time_series = running_mean(core_time_series,N)
        ax2.plot(uniq_time_sorted,tmp_core_time_series,lw=lws2[dum],ls='solid',color=colors2[dum])
        dum+=1    

    # theta_pert    
    dum = 0
    for sim in sims3:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        core_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            core_time_series.append(np.size(tmpid[0]))
        core_time_series = np.array(core_time_series)
        tmp_core_time_series = running_mean(core_time_series,N)
        ax3.plot(uniq_time_sorted,tmp_core_time_series,lw=lws3[dum],ls=lss3[dum],color=colors3[dum])
        dum+=1  


    #===============================
    # Mean Core Area
    #===============================
    # STOCH_SHORT
    dum = 0
    for sim in sims1:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_area = area[sim]
        tmp_area = np.array(tmp_area)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        mean_area_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_mean_area = np.nanmean(tmp_area[tmpid[0]])
            mean_area_time_series.append(tmp_mean_area)
        mean_area_time_series = np.array(mean_area_time_series)
        print(aaaa)
        tmp_mean_area_time_series = running_mean(mean_area_time_series,N)
        ax4.plot(uniq_time_sorted,tmp_mean_area_time_series,lw=lws1[dum],ls='solid',color=colors1[dum])
        dum+=1

    # STOCH_LONG    
    dum = 0
    for sim in sims2:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_area = area[sim]
        tmp_area = np.array(tmp_area)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        mean_area_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_mean_area = np.nanmean(tmp_area[tmpid[0]])
            mean_area_time_series.append(tmp_mean_area)
        mean_area_time_series = np.array(mean_area_time_series)
        tmp_mean_area_time_series = running_mean(mean_area_time_series,N)
        ax5.plot(uniq_time_sorted,tmp_mean_area_time_series,lw=lws2[dum],ls='solid',color=colors2[dum])
        dum+=1 

    # theta_pert    
    dum = 0
    for sim in sims3:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_area = area[sim]
        tmp_area = np.array(tmp_area)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        mean_area_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_mean_area = np.nanmean(tmp_area[tmpid[0]])
            mean_area_time_series.append(tmp_mean_area)
        mean_area_time_series = np.array(mean_area_time_series)
        tmp_mean_area_time_series = running_mean(mean_area_time_series,N)
        ax6.plot(uniq_time_sorted,tmp_mean_area_time_series,lw=lws3[dum],ls=lss3[dum],color=colors3[dum])
        dum+=1


    #===============================
    # Total Core Area
    #===============================
    # STOCH_SHORT
    dum = 0
    for sim in sims1:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_area = area[sim]
        tmp_area = np.array(tmp_area)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        tot_area_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_tot_area = np.nansum(tmp_area[tmpid[0]])
            tot_area_time_series.append(tmp_tot_area)
        tot_area_time_series = np.array(tot_area_time_series)
        tmp_tot_area_time_series = running_mean(tot_area_time_series,N)
        ax7.plot(uniq_time_sorted,tmp_tot_area_time_series,lw=lws1[dum],ls='solid',color=colors1[dum])
        dum+=1

    # STOCH_LONG    
    dum = 0
    for sim in sims2:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_area = area[sim]
        tmp_area = np.array(tmp_area)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        tot_area_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_to_area = np.nansum(tmp_area[tmpid[0]])
            tot_area_time_series.append(tmp_to_area)
        tot_area_time_series = np.array(tot_area_time_series)
        tmp_tot_area_time_series = running_mean(tot_area_time_series,N)
        ax8.plot(uniq_time_sorted,tmp_tot_area_time_series,lw=lws2[dum],ls='solid',color=colors2[dum])
        dum+=1 

    # theta_pert    
    dum = 0
    for sim in sims3:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_area = area[sim]
        tmp_area = np.array(tmp_area)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        tot_area_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_tot_area = np.nansum(tmp_area[tmpid[0]])
            tot_area_time_series.append(tmp_tot_area)
        tot_area_time_series = np.array(tot_area_time_series)
        tmp_tot_area_time_series = running_mean(tot_area_time_series,N)
        ax9.plot(uniq_time_sorted,tmp_tot_area_time_series,lw=lws3[dum],ls=lss3[dum],color=colors3[dum])
        dum+=1    

    #=======================================
    # Domain Median Multiplicative Factor
    #=======================================
    tmp_time = sorted(np.unique(time['baseline']))


    if False:
        # STOCH_SHORT
        dum = 0
        for sim in sims4:
            ax10.plot(tmp_time,median_mult_fac_dict[sim],lw=2,ls='solid',color=colors4[dum])
            dum+=1

        # STOCH_LONG    
        dum = 0
        for sim in sims5:
            ax11.plot(tmp_time,median_mult_fac_dict[sim],lw=2,ls='solid',color=colors5[dum])
            dum+=1    





    axlist1 = [ax1,ax2,ax3]
    for ax in axlist1:
        ax.grid(which='both',ls='solid',lw=1,color='grey')
        ax.tick_params(labelsize=23)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Number of CCR \nCores',fontsize=23)
        ax.set_xlabel('Time [UTC]',fontsize=23)
        ax.set_xticks([tmp_time[0],tmp_time[24],tmp_time[48],\
                            tmp_time[72],tmp_time[96]])
        ax.set_ylim(0,100)

    axlist2 = [ax4,ax5,ax6]
    for ax in axlist2:
        ax.grid(which='both',ls='solid',lw=1,color='grey')
        ax.tick_params(labelsize=23)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Mean CCR Core\nArea [km$^{2}$]',fontsize=23)
        ax.set_xlabel('Time [UTC]',fontsize=23)
        ax.set_xticks([tmp_time[0],tmp_time[24],tmp_time[48],\
                            tmp_time[72],tmp_time[96]])
        ax.set_ylim(0,225)


    axlist3 = [ax7,ax8,ax9]
    for ax in axlist3:
        ax.grid(which='both',ls='solid',lw=1,color='grey')
        ax.tick_params(labelsize=23)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Total CCR Core \nArea [km$^{2}$]',fontsize=23)
        ax.set_xlabel('Time [UTC]',fontsize=23)
        ax.set_xticks([tmp_time[0],tmp_time[24],tmp_time[48],\
                            tmp_time[72],tmp_time[96]])
        ax.set_ylim(0,8000)    

    #axlist4 = [ax10,ax11]
    #for ax in axlist4:
    if False:
        ax.grid(which='both',ls='solid',lw=1,color='grey')
        ax.tick_params(labelsize=20)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Domain Median $F$',fontsize=20)
        ax.set_xlabel('Time [UTC]',fontsize=20)
        ax.set_xticks([tmp_time[0],tmp_time[24],tmp_time[48],\
                            tmp_time[72],tmp_time[96]])
        ax.axhline(1,lw=4,color='black')
        ax.set_yscale('log')
        ax.set_ylim(0.25,4)
        #ax.set_ylim(08000)  


    plt.subplots_adjust(hspace=0.4,wspace=0.41)    


    plt.text(0.5,1.1,'STOCH\_SHORT',fontsize=33,fontweight='bold',\
                bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
                ha='center')
    plt.text(0.5,1.1,'STOCH\_LONG',fontsize=33,fontweight='bold',\
                 bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
                 ha='center')
    plt.text(0.5,1.1,'$\\theta$-pert',fontsize=33,fontweight='bold',\
                 bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
                 ha='center')


    #axlist = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11]
    axlist = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]
    #labs = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)']
    labs = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)']
    dum = 0
    for ax in axlist:
        ax.text(0.05,0.77,labs[dum],fontsize=50,fontweight='bold',transform=ax.transAxes)
        dum+=1

    custom_lines = [Line2D([0],[0],color = 'black',lw=6),\
                        Line2D([0],[0],color='grey',lw=6),\
                        Line2D([0],[0],color='black',ls='dotted',lw=6),\
                        Line2D([0],[0],color='black',ls='dashed',lw=6),\
                        Line2D([0],[0],color='darkorange',ls='-.',lw=6),\
                       ]

    ax8.legend(custom_lines,['BASELINE','SPOL','NO\_MIXING','4X','BASELINE\_1KM\_CG'],fontsize=30,\
                        loc='lower center',bbox_to_anchor=(0.5,-0.9),ncol=3)

    

    #print(aaaaa)    
    








if True:
    #----------------------------------
    #----------------------------------
    #----------------------------------
    #----------------------------------
    # Time Series of Core Number
    #----------------------------------
    #----------------------------------
    #----------------------------------
    #----------------------------------
    import matplotlib.dates as mdates
    plt.rc('text',usetex=True)
    dfmt = mdates.DateFormatter('%d-%H')

    sims1 = ['baseline',\
             'stoch1_short',\
             'stoch2_short',\
             'stoch3_short',\
             'stoch4_short',\
             'stoch5_short',\
             'SPOL',\
            ]

    sims4 = ['stoch1_short',\
             'stoch2_short',\
             'stoch3_short',\
             'stoch4_short',\
             'stoch5_short']

    sims2 = ['baseline',\
             'stoch1_long',\
             'stoch2_long',\
             'stoch3_long',\
             'stoch4_long',\
             'stoch5_long',\
             'SPOL']

    sims5 = ['stoch1_long',\
             'stoch2_long',\
             'stoch3_long',\
             'stoch4_long',\
             'stoch5_long']         

    sims3 = ['baseline',\
             'theta_pert1',\
             'theta_pert2',\
             'theta_pert3',\
             'theta_pert4',\
             'theta_pert5',\
             'SPOL',\
             'no_mixing',\
             '4x',\
             'baseline_1km_CG',\
            ]

    colors1 = ['black',\
             'lightsalmon',\
             'salmon',\
             'red',\
             'darkred',\
             'firebrick',\
              'grey']

    colors4 = ['lightsalmon',\
             'salmon',\
             'red',\
             'darkred',\
             'firebrick']

    colors2 = ['black',\
             'powderblue',\
             'deepskyblue',\
             'dodgerblue',\
             'blue',\
             'navy',\
             'grey',\
             ]  

    colors5 = ['powderblue',\
             'deepskyblue',\
             'dodgerblue',\
             'blue',\
             'navy'] 

    colors3 = ['black',\
             'orchid',\
             'fuchsia',\
             'mediumorchid',\
             'purple',\
             'darkorchid',\
             'grey',\
             'black',\
             'black',\
             'darkorange',\
             ]  

    lws1 = [5,3,3,3,3,3,5]
    lws2 = [5,3,3,3,3,3,5]
    lws3 = [5,3,3,3,3,3,5,5,5,5]
    lss3 = ['solid','solid','solid','solid','solid','solid','solid','dotted','dashed','-.']

    fig = plt.figure(figsize=(22,16))

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    #axlist = [ax1,ax2,ax3,ax4,ax5,ax6]

    #===============================
    # Number of Cores
    #===============================
    # STOCH_SHORT
    dum = 0
    for sim in sims1:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        core_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            core_time_series.append(np.size(tmpid[0]))
        core_time_series = np.array(core_time_series)
        tmp_core_time_series = running_mean(core_time_series,N)
        ax1.plot(uniq_time_sorted,tmp_core_time_series,lw=lws1[dum],ls='solid',color=colors1[dum])
        dum+=1

    # STOCH_LONG    
    dum = 0
    for sim in sims2:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        core_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            core_time_series.append(np.size(tmpid[0]))
        core_time_series = np.array(core_time_series)
        tmp_core_time_series = running_mean(core_time_series,N)
        ax2.plot(uniq_time_sorted,tmp_core_time_series,lw=lws2[dum],ls='solid',color=colors2[dum])
        dum+=1    


    #===============================
    # Mean Core Area
    #===============================
    # STOCH_SHORT
    dum = 0
    for sim in sims1:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_area = area[sim]
        tmp_area = np.array(tmp_area)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        mean_area_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_mean_area = np.nanmean(tmp_area[tmpid[0]])
            mean_area_time_series.append(tmp_mean_area)
        mean_area_time_series = np.array(mean_area_time_series)
        tmp_mean_area_time_series = running_mean(mean_area_time_series,N)
        ax3.plot(uniq_time_sorted,tmp_mean_area_time_series,lw=lws1[dum],ls='solid',color=colors1[dum])
        dum+=1

    # STOCH_LONG    
    dum = 0
    for sim in sims2:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_area = area[sim]
        tmp_area = np.array(tmp_area)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        mean_area_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_mean_area = np.nanmean(tmp_area[tmpid[0]])
            mean_area_time_series.append(tmp_mean_area)
        mean_area_time_series = np.array(mean_area_time_series)
        tmp_mean_area_time_series = running_mean(mean_area_time_series,N)
        ax4.plot(uniq_time_sorted,tmp_mean_area_time_series,lw=lws2[dum],ls='solid',color=colors2[dum])
        dum+=1 

    tmp_time = uniq_time_sorted

    axlist1 = [ax1,ax2]
    for ax in axlist1:
        ax.grid(which='both',ls='solid',lw=2,color='grey')
        ax.tick_params(labelsize=45)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Number of CCRs',fontsize=45)
        ax.set_xlabel('Time [UTC]',fontsize=45)
        ax.set_xticks([tmp_time[0],tmp_time[24],tmp_time[48],\
                            tmp_time[72],tmp_time[96]])
        ax.set_ylim(0,100)

    axlist2 = [ax3,ax4]
    for ax in axlist2:
        ax.grid(which='both',ls='solid',lw=2,color='grey')
        ax.tick_params(labelsize=45)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Mean CCR Area [km$^{2}$]',fontsize=45)
        ax.set_xlabel('Time [UTC]',fontsize=45)
        ax.set_xticks([tmp_time[0],tmp_time[24],tmp_time[48],\
                            tmp_time[72],tmp_time[96]])
        ax.set_ylim(0,225)


    plt.subplots_adjust(hspace=0.4,wspace=0.3)    


    plt.text(0.5,1.1,'STOCH\_SHORT',fontsize=60,fontweight='bold',\
                bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
                ha='center')
    plt.text(0.5,1.1,'STOCH\_LONG',fontsize=60,fontweight='bold',\
                 bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
                 ha='center')




    custom_lines = [Line2D([0],[0],color = 'black',lw=6),\
                        Line2D([0],[0],color='grey',lw=6),\
                       ]

    ax3.legend(custom_lines,['BASELINE','SPOL'],fontsize=50,\
                        loc='lower center',bbox_to_anchor=(1,-0.75),ncol=2)


    #print(aaaaa)    
    


#----------------------------------
#----------------------------------
#----------------------------------
#----------------------------------
# Time Series of Core-averaged
# properties, particularly ETH
# and rain rate
#----------------------------------
#----------------------------------
#----------------------------------
#----------------------------------
import matplotlib.dates as mdates
plt.rc('text',usetex=True)
dfmt = mdates.DateFormatter('%d-%H')


if False:

    sims1 = ['baseline',\
             'stoch1_short',\
             'stoch2_short',\
             'stoch3_short',\
             'stoch4_short',\
             'stoch5_short',\
             'SPOL',\
            ]


    sims2 = ['baseline',\
             'stoch1_long',\
             'stoch2_long',\
             'stoch3_long',\
             'stoch4_long',\
             'stoch5_long',\
             'SPOL']

    sims3 = ['baseline',\
             'theta_pert1',\
             'theta_pert2',\
             'theta_pert3',\
             'theta_pert4',\
             'theta_pert5',\
             'SPOL',\
             'no_mixing',\
             '4x',\
            ]

    colors1 = ['black',\
             'lightsalmon',\
             'salmon',\
             'red',\
             'darkred',\
             'firebrick',\
              'grey']


    colors2 = ['black',\
             'powderblue',\
             'deepskyblue',\
             'dodgerblue',\
             'blue',\
             'navy',\
             'grey',\
             ]  


    colors3 = ['black',\
             'orchid',\
             'fuchsia',\
             'mediumorchid',\
             'purple',\
             'darkorchid',\
             'grey',\
             'black',\
             'black',\
             ]  

    lws1 = [4,2,2,2,2,2,4]
    lws2 = [4,2,2,2,2,2,4]
    lws3 = [4,2,2,2,2,2,4,4,4]
    lss3 = ['solid','solid','solid','solid','solid','solid','solid','dotted','dashed']

    fig = plt.figure(figsize=(18,16))

    ax1 = fig.add_subplot(431)
    ax2 = fig.add_subplot(432)
    ax3 = fig.add_subplot(433)
    ax4 = fig.add_subplot(434)
    ax5 = fig.add_subplot(435)
    ax6 = fig.add_subplot(436)
    ax7 = fig.add_subplot(437)
    ax8 = fig.add_subplot(438)
    ax9 = fig.add_subplot(439)
    ax10 = fig.add_subplot(4,3,10)
    ax11 = fig.add_subplot(4,3,11)

    #axlist = [ax1,ax2,ax3,ax4,ax5,ax6]




    #===============================
    # Mean Core Average Rain Rate
    #===============================

    # STOCH_SHORT
    dum = 0
    for sim in sims1:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_avg_rain_rate = avg_rain_rate[sim]
        tmp_avg_rain_rate = np.array(tmp_avg_rain_rate)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        mean_avg_rain_rate_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_mean_avg_rain_rate = np.nanmean(tmp_avg_rain_rate[tmpid[0]])
            mean_avg_rain_rate_time_series.append(tmp_mean_avg_rain_rate)
        mean_avg_rain_rate_time_series = np.array(mean_avg_rain_rate_time_series)
        tmp_mean_avg_rain_rate_time_series = running_mean(mean_avg_rain_rate_time_series,N)
        ax1.plot(uniq_time_sorted,tmp_mean_avg_rain_rate_time_series,lw=lws1[dum],ls='solid',color=colors1[dum])
        dum+=1

    # STOCH_LONG    
    dum = 0
    for sim in sims2:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_avg_rain_rate = avg_rain_rate[sim]
        tmp_avg_rain_rate = np.array(tmp_avg_rain_rate)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        mean_avg_rain_rate_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_mean_avg_rain_rate = np.nanmean(tmp_avg_rain_rate[tmpid[0]])
            mean_avg_rain_rate_time_series.append(tmp_mean_avg_rain_rate)
        mean_avg_rain_rate_time_series = np.array(mean_avg_rain_rate_time_series)
        tmp_mean_avg_rain_rate_time_series = running_mean(mean_avg_rain_rate_time_series,N)
        ax2.plot(uniq_time_sorted,tmp_mean_avg_rain_rate_time_series,lw=lws2[dum],ls='solid',color=colors2[dum])
        dum+=1 

    # theta_pert    
    dum = 0
    for sim in sims3:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_avg_rain_rate = avg_rain_rate[sim]
        tmp_avg_rain_rate = np.array(tmp_avg_rain_rate)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        mean_avg_rain_rate_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_mean_avg_rain_rate = np.nanmean(tmp_avg_rain_rate[tmpid[0]])
            mean_avg_rain_rate_time_series.append(tmp_mean_avg_rain_rate)
        mean_avg_rain_rate_time_series = np.array(mean_avg_rain_rate_time_series)
        tmp_mean_avg_rain_rate_time_series = running_mean(mean_avg_rain_rate_time_series,N)
        ax3.plot(uniq_time_sorted,tmp_mean_avg_rain_rate_time_series,lw=lws3[dum],ls=lss3[dum],color=colors3[dum])
        dum+=1

    print(aaaaaa)
    #===============================
    # Total Core Area
    #===============================
    # STOCH_SHORT
    dum = 0
    for sim in sims1:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_area = area[sim]
        tmp_area = np.array(tmp_area)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        tot_area_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_tot_area = np.nansum(tmp_area[tmpid[0]])
            tot_area_time_series.append(tmp_tot_area)
        tot_area_time_series = np.array(tot_area_time_series)
        tmp_tot_area_time_series = running_mean(tot_area_time_series,N)
        ax7.plot(uniq_time_sorted,tmp_tot_area_time_series,lw=lws1[dum],ls='solid',color=colors1[dum])
        dum+=1

    # STOCH_LONG    
    dum = 0
    for sim in sims2:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_area = area[sim]
        tmp_area = np.array(tmp_area)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        tot_area_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_to_area = np.nansum(tmp_area[tmpid[0]])
            tot_area_time_series.append(tmp_to_area)
        tot_area_time_series = np.array(tot_area_time_series)
        tmp_tot_area_time_series = running_mean(tot_area_time_series,N)
        ax8.plot(uniq_time_sorted,tmp_tot_area_time_series,lw=lws2[dum],ls='solid',color=colors2[dum])
        dum+=1 

    # theta_pert    
    dum = 0
    for sim in sims3:
        tmp_time = convective_objects['time'][sim]
        tmp_time = np.array(tmp_time)
        tmp_area = area[sim]
        tmp_area = np.array(tmp_area)

        uniq_time = np.unique(tmp_time)
        uniq_time_sorted = sorted(uniq_time)
        nt = len(uniq_time_sorted)
        tot_area_time_series = []
        for tt in range(nt):
            tmpid = np.where(tmp_time == uniq_time_sorted[tt])
            tmp_tot_area = np.nansum(tmp_area[tmpid[0]])
            tot_area_time_series.append(tmp_tot_area)
        tot_area_time_series = np.array(tot_area_time_series)
        tmp_tot_area_time_series = running_mean(tot_area_time_series,N)
        ax9.plot(uniq_time_sorted,tmp_tot_area_time_series,lw=lws3[dum],ls=lss3[dum],color=colors3[dum])
        dum+=1    

    #=======================================
    # Domain Median Multiplicative Factor
    #=======================================

    tmp_time = sorted(np.unique(time['baseline']))
    # STOCH_SHORT
    dum = 0
    for sim in sims4:
        ax10.plot(tmp_time,median_mult_fac_dict[sim],lw=2,ls='solid',color=colors4[dum])
        dum+=1

    # STOCH_LONG    
    dum = 0
    for sim in sims5:
        ax11.plot(tmp_time,median_mult_fac_dict[sim],lw=2,ls='solid',color=colors5[dum])
        dum+=1    





    axlist1 = [ax1,ax2,ax3]
    for ax in axlist1:
        ax.grid(which='both',ls='solid',lw=1,color='grey')
        ax.tick_params(labelsize=20)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Number of Cores',fontsize=20)
        ax.set_xlabel('Time [UTC]',fontsize=20)
        ax.set_xticks([tmp_time[0],tmp_time[24],tmp_time[48],\
                            tmp_time[72],tmp_time[96]])
        ax.set_ylim(0,80)

    axlist2 = [ax4,ax5,ax6]
    for ax in axlist2:
        ax.grid(which='both',ls='solid',lw=1,color='grey')
        ax.tick_params(labelsize=20)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Mean Core Area [km$^{2}$]',fontsize=20)
        ax.set_xlabel('Time [UTC]',fontsize=20)
        ax.set_xticks([tmp_time[0],tmp_time[24],tmp_time[48],\
                            tmp_time[72],tmp_time[96]])
        ax.set_ylim(0,225)


    axlist3 = [ax7,ax8,ax9]
    for ax in axlist3:
        ax.grid(which='both',ls='solid',lw=1,color='grey')
        ax.tick_params(labelsize=20)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Total Core Area [km$^{2}$]',fontsize=20)
        ax.set_xlabel('Time [UTC]',fontsize=20)
        ax.set_xticks([tmp_time[0],tmp_time[24],tmp_time[48],\
                            tmp_time[72],tmp_time[96]])
        ax.set_ylim(0,8000)    

    axlist4 = [ax10,ax11]
    for ax in axlist4:
        ax.grid(which='both',ls='solid',lw=1,color='grey')
        ax.tick_params(labelsize=20)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Domain Median $F$',fontsize=20)
        ax.set_xlabel('Time [UTC]',fontsize=20)
        ax.set_xticks([tmp_time[0],tmp_time[24],tmp_time[48],\
                            tmp_time[72],tmp_time[96]])
        ax.axhline(1,lw=4,color='black')
        ax.set_yscale('log')
        ax.set_ylim(0.25,4)
        #ax.set_ylim(08000)  


    plt.subplots_adjust(hspace=0.4,wspace=0.3)    


    plt.text(0.5,1.1,'STOCH\_SHORT',fontsize=33,fontweight='bold',\
                bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
                ha='center')
    plt.text(0.5,1.1,'STOCH\_LONG',fontsize=33,fontweight='bold',\
                 bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
                 ha='center')
    plt.text(0.5,1.1,'$\\theta$-pert',fontsize=33,fontweight='bold',\
                 bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
                 ha='center')


    axlist = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11]
    labs = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)']
    dum = 0
    for ax in axlist:
        ax.text(0.05,0.77,labs[dum],fontsize=50,fontweight='bold',transform=ax.transAxes)
        dum+=1

    custom_lines = [Line2D([0],[0],color = 'black',lw=6),\
                        Line2D([0],[0],color='grey',lw=6),\
                        Line2D([0],[0],color='black',ls='dotted',lw=6),\
                        Line2D([0],[0],color='black',ls='dashed',lw=6),\
                       ]

    ax9.legend(custom_lines,['BASELINE','SPOL','NO\_MIXING','4X'],fontsize=30,\
                        loc='lower center',bbox_to_anchor=(0.5,-1.6),ncol=1)


    print(aaaaa)    

















#----------------------------------
# PDF of cell area
#----------------------------------

area_pdf_run = 0
if area_pdf_run == 1:

    fig = plt.figure(figsize=(18,10))

    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
    for ax in axlist:
        ax.grid(which='both',ls='solid',color='grey',lw=0.5)
        ax.set_xlabel('Cell Area [km$^{2}$]',fontsize=15)
        ax.tick_params(labelsize=15)
        ax.set_xscale('log')
        ax.set_xlim(20,0.5e3)
        ax.set_ylabel('Frequency',fontsize=15)
        #ax.set_ylabel('Count',fontsize=15)

    
    bins = np.arange(0,27*300,18)
    midbins = []
    for ii in range(len(bins)-1):
        midbins.append(0.5*(bins[ii]+bins[ii+1]))
    midbins = np.array(midbins)

    keys1 = ['baseline','stoch1_short','stoch2_short','stoch3_short',\
            'stoch4_short','stoch5_short','SPOL','stoch1_short_above500m_sigma1',\
            'stoch2_short_above500m_sigma1','stoch1_short_mynn','baseline_1km_CG']
    
    keys2 = ['baseline','stoch1_long','stoch2_long','stoch3_long',\
            'stoch4_long','stoch5_long','SPOL','stoch1_short_mynn','baseline_1km_CG']   
    
    keys3 = ['baseline','theta_pert1','theta_pert2','theta_pert3',\
            'theta_pert4','theta_pert5','SPOL','stoch1_short_mynn','baseline_1km_CG']      

    
    colors1 = ['black',\
               'lightsalmon',\
               'salmon',\
               'red',\
               'darkred',\
               'firebrick','grey','green','lime','gold','darkorange']
    
    colors2 = ['black',\
         'powderblue',\
         'deepskyblue',\
         'dodgerblue',\
         'blue',\
         'navy','grey','gold','darkorange']
    
    colors3 = ['black',\
             'orchid',\
             'fuchsia',\
             'mediumorchid',\
             'purple',\
             'darkorchid','grey','gold','darkorange']
    

    dum = 0
    for key in keys1:
        val = area[key]
        val = np.array(val)
        val_10dbz = np.array(max_echo_top_10[key])
        if key == 'SPOL':
            val_10dbz = val_10dbz*1000.
        #tmpid = np.where(val_10dbz >= 8000.)
        #tmpid = np.where((val_10dbz >= 6000.) & (val_10dbz < 8000.))
        
        #val = val[tmpid]
        val = val[val > 9.]
        
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax1.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
        #tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax1.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=2,ls='solid')#,marker=markers[dum])
        #ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=2,ls='solid')#,marker=markers[dum])#
        dum = dum + 1    
        
    dum = 0
    for key in keys2:
        val = area[key]
        val = np.array(val) 
        val_10dbz = np.array(max_echo_top_10[key])
        if key == 'SPOL':
            val_10dbz = val_10dbz*1000.       
        #tmpid = np.where(val_10dbz >= 8000.)
        #tmpid = np.where((val_10dbz >= 6000.) & (val_10dbz < 8000.))
        #val = val[tmpid]        
        val = val[val > 9.]
        
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax2.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax5.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
        #tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax2.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
        #         lw=2,ls='solid')#,marker=markers[dum])
        #ax5.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
        #         lw=2,ls='solid')#,marker=markers[dum])#
        dum = dum + 1    
        
    dum = 0
    for key in keys3:
        val = area[key]
        val = np.array(val)
        val_10dbz = np.array(max_echo_top_10[key])
        #print(np.max(val_10dbz))
        if key == 'SPOL':
            val_10dbz = val_10dbz*1000.
        #tmpid = np.where(val_10dbz >= 8000.)
        #tmpid = np.where((val_10dbz >= 6000.) & (val_10dbz < 8000.))
        #val = val[tmpid]              
        val = val[val > 9.]
        
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax6.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
        #tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
        #         lw=2,ls='solid')#,marker=markers[dum])
        #ax6.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
        #         lw=2,ls='solid')#,marker=markers[dum])#
        dum = dum + 1            
   

    #ax3.legend(fontsize=15,bbox_to_anchor=(2.,-0.18),ncol=4)
    ax4.legend(fontsize=15,bbox_to_anchor=(0.6,-1.5),ncol=1,loc='lower center')
    ax5.legend(fontsize=15,bbox_to_anchor=(0.6,-1.3),ncol=1,loc='lower center')
    ax6.legend(fontsize=15,bbox_to_anchor=(0.6,-1.3),ncol=1,loc='lower center')

    ax4.set_yscale('log')
    ax5.set_yscale('log')
    ax6.set_yscale('log')
    #ax2.set_yscale('log')
    #ax4.set_yscale('log')
    #ax6.set_yscale('log')
    #ax8.set_yscale('log')
    #ax1.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax2.set_title('PDF, $Log-y$',fontsize=25)
    #ax3.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax4.set_title('Histogram, $Log-y$',fontsize=25)
    #ax5.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax6.set_title('PDF, $Log-y$',fontsize=25)
    #ax7.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax8.set_title('Histogram, $Log-y$',fontsize=25)
    
    #plt.suptitle('8 Dec.',fontsize=40)
    #plt.suptitle('All Times, 6 km <= 10 dBZ ETH < 8 km',fontsize=40)
    #plt.suptitle('Before 21 UTC on 7th',fontsize=40)
    plt.suptitle('All Times',fontsize=40)
    
    
    plt.subplots_adjust(wspace=0.3,hspace=0.35)
    
    
    plt.show()
    print(aaaaaa)

    
#----------------------------------
# PDF of cell max 0 dBZ ETH
#----------------------------------

cellmax_0dbz_run = 0
if cellmax_0dbz_run == 1:

    fig = plt.figure(figsize=(18,10))

    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
    for ax in axlist:
        ax.grid(which='both',ls='solid',color='grey',lw=0.5)
        ax.set_xlabel('Core Max. 0 dBZ ETH [km]',fontsize=15)
        ax.tick_params(labelsize=15)
        #ax.set_xscale('log')
        #ax.set_xlim(20,0.5e3)
        #ax.set_ylabel('Frequency',fontsize=15)
        ax.set_ylabel('Frequency',fontsize=15)

    
    #bins = np.arange(0,21,1)
    bins = np.arange(0,21,1)
    midbins = []
    for ii in range(len(bins)-1):
        midbins.append(0.5*(bins[ii]+bins[ii+1]))
    midbins = np.array(midbins)

    keys1 = ['baseline','stoch1_short','stoch2_short','stoch3_short',\
            'stoch4_short','stoch5_short','SPOL','stoch1_short_above500m_sigma1',\
            'stoch2_short_above500m_sigma1','stoch1_short_mynn','baseline_1km_CG']
    keys2 = ['baseline','stoch1_long','stoch2_long','stoch3_long',\
            'stoch4_long','stoch5_long','SPOL','stoch1_short_mynn','baseline_1km_CG']   
    keys3 = ['baseline','theta_pert1','theta_pert2','theta_pert3',\
            'theta_pert4','theta_pert5','SPOL','stoch1_short_mynn','baseline_1km_CG']      

    
    colors1 = ['black',\
               'lightsalmon',\
               'salmon',\
               'red',\
               'darkred',\
               'firebrick','grey','green',\
              'lime','gold','darkorange']
    
    colors2 = ['black',\
         'powderblue',\
         'deepskyblue',\
         'dodgerblue',\
         'blue',\
         'navy','grey','gold','darkorange']
    
    colors3 = ['black',\
             'orchid',\
             'fuchsia',\
             'mediumorchid',\
             'purple',\
             'darkorchid','grey','gold','darkorange']
    

    dum = 0
    for key in keys1:
        val = max_echo_top_0[key]
        val = np.array(val)
        if key != 'SPOL':
            val = val/1.e3
        # filter
        #val = val[val > 6.]
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax1.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
       # tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])
        #ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])#
        dum = dum + 1    
        
    dum = 0
    for key in keys2:
        val = max_echo_top_0[key]
        val = np.array(val)     
        if key != 'SPOL':
            val = val/1.e3
        # filter
        #val = val[val > 6.]            
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax2.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax5.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
       # tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])
        #ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])#
        dum = dum + 1    
        
    dum = 0
    for key in keys3:
        val = max_echo_top_0[key]
        val = np.array(val)
        if key != 'SPOL':
            val = val/1.e3
        # filter
       # val = val[val > 6.]            
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax6.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
       # tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])
        #ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])#
        dum = dum + 1            
   

    #ax3.legend(fontsize=15,bbox_to_anchor=(2.,-0.18),ncol=4)
    ax4.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)
    ax5.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)
    ax6.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)

    ax4.set_yscale('log')
    ax5.set_yscale('log')
    ax6.set_yscale('log')
    #ax2.set_yscale('log')
    #ax4.set_yscale('log')
    #ax6.set_yscale('log')
    #ax8.set_yscale('log')
    #ax1.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax2.set_title('PDF, $Log-y$',fontsize=25)
    #ax3.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax4.set_title('Histogram, $Log-y$',fontsize=25)
    #ax5.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax6.set_title('PDF, $Log-y$',fontsize=25)
    #ax7.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax8.set_title('Histogram, $Log-y$',fontsize=25)
    
    
    plt.subplots_adjust(wspace=0.3,hspace=0.35)
    
    plt.suptitle('Before 21 UTC on 7th',fontsize=40)
    #plt.suptitle('All Times',fontsize=40)
    #plt.suptitle('8 Dec.',fontsize=40)
    
    plt.show()
    print(aaaaaa)

    
    
#----------------------------------
# Histogram of cell max 0 dBZ ETH
#----------------------------------

cellmax_0dbz_run_hist = 0
if cellmax_0dbz_run_hist == 1:

    fig = plt.figure(figsize=(18,10))

    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
    for ax in axlist:
        ax.grid(which='both',ls='solid',color='grey',lw=0.5)
        ax.set_xlabel('Core Max. 0 dBZ ETH [km]',fontsize=15)
        ax.tick_params(labelsize=15)
        #ax.set_xscale('log')
        #ax.set_xlim(20,0.5e3)
        #ax.set_ylabel('Frequency',fontsize=15)
        ax.set_ylabel('Count',fontsize=15)

    
    #bins = np.arange(0,21,1)
    bins = np.arange(0,21,1)
    midbins = []
    for ii in range(len(bins)-1):
        midbins.append(0.5*(bins[ii]+bins[ii+1]))
    midbins = np.array(midbins)

    keys1 = ['baseline','stoch1_short','stoch2_short','stoch3_short',\
            'stoch4_short','stoch5_short','SPOL','stoch1_short_above500m_sigma1',\
            'stoch2_short_above500m_sigma1','stoch1_short_mynn','baseline_1km_CG']
    keys2 = ['baseline','stoch1_long','stoch2_long','stoch3_long',\
            'stoch4_long','stoch5_long','SPOL','stoch1_short_mynn','baseline_1km_CG']   
    keys3 = ['baseline','theta_pert1','theta_pert2','theta_pert3',\
            'theta_pert4','theta_pert5','SPOL','stoch1_short_mynn','baseline_1km_CG']      

    
    colors1 = ['black',\
               'lightsalmon',\
               'salmon',\
               'red',\
               'darkred',\
               'firebrick','grey','green',\
              'lime','gold','darkorange']
    
    colors2 = ['black',\
         'powderblue',\
         'deepskyblue',\
         'dodgerblue',\
         'blue',\
         'navy','grey','gold','darkorange']
    
    colors3 = ['black',\
             'orchid',\
             'fuchsia',\
             'mediumorchid',\
             'purple',\
             'darkorchid','grey','gold','darkorange']
    

    dum = 0
    for key in keys1:
        val = max_echo_top_0[key]
        val = np.array(val)
        if key != 'SPOL':
            val = val/1.e3
        # filter
        #val = val[val > 6.]

        # Histogram
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax1.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
                 lw=2)#,marker=markers[dum])
        ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
                 lw=2)#,marker=markers[dum])#
        dum = dum + 1    
        
    dum = 0
    for key in keys2:
        val = max_echo_top_0[key]
        val = np.array(val)     
        if key != 'SPOL':
            val = val/1.e3
        # filter
        #val = val[val > 6.]            

        # Histogram
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax2.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
                 lw=2)#,marker=markers[dum])
        ax5.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
                 lw=2)#,marker=markers[dum])#
        dum = dum + 1    
        
    dum = 0
    for key in keys3:
        val = max_echo_top_0[key]
        val = np.array(val)
        if key != 'SPOL':
            val = val/1.e3
        # filter
       # val = val[val > 6.]            

        # Histogram
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
                 lw=2)#,marker=markers[dum])
        ax6.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
                 lw=2)#,marker=markers[dum])#
        dum = dum + 1            
   

    #ax3.legend(fontsize=15,bbox_to_anchor=(2.,-0.18),ncol=4)
    ax4.legend(fontsize=15,bbox_to_anchor=(0.6,-1.5),ncol=1,loc='lower center')
    ax5.legend(fontsize=15,bbox_to_anchor=(0.6,-1.3),ncol=1,loc='lower center')
    ax6.legend(fontsize=15,bbox_to_anchor=(0.6,-1.3),ncol=1,loc='lower center')

    ax4.set_yscale('log')
    ax5.set_yscale('log')
    ax6.set_yscale('log')
    #ax2.set_yscale('log')
    #ax4.set_yscale('log')
    #ax6.set_yscale('log')
    #ax8.set_yscale('log')
    #ax1.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax2.set_title('PDF, $Log-y$',fontsize=25)
    #ax3.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax4.set_title('Histogram, $Log-y$',fontsize=25)
    #ax5.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax6.set_title('PDF, $Log-y$',fontsize=25)
    #ax7.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax8.set_title('Histogram, $Log-y$',fontsize=25)
    
    
    plt.subplots_adjust(wspace=0.3,hspace=0.35)
    
    #plt.suptitle('Before 21 UTC on 7th',fontsize=40)
    plt.suptitle('All Times',fontsize=40)
    #plt.suptitle('8 Dec.',fontsize=40)
    
    plt.show()
    print(aaaaaa)


#----------------------------------
# PDF of cell max 20 dBZ ETH
#----------------------------------

cellmax_20dbz_run = 0
if cellmax_20dbz_run == 1:

    fig = plt.figure(figsize=(18,10))

    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
    for ax in axlist:
        ax.grid(which='both',ls='solid',color='grey',lw=0.5)
        ax.set_xlabel('Core Max. 20 dBZ ETH [km]',fontsize=15)
        ax.tick_params(labelsize=15)
        #ax.set_xscale('log')
        #ax.set_xlim(20,0.5e3)
        ax.set_ylabel('Frequency',fontsize=15)

    
    #bins = np.arange(0,21,1)
    bins = np.arange(0,21,1)
    midbins = []
    for ii in range(len(bins)-1):
        midbins.append(0.5*(bins[ii]+bins[ii+1]))
    midbins = np.array(midbins)

    keys1 = ['baseline','stoch1_short','stoch2_short','stoch3_short',\
            'stoch4_short','stoch5_short','SPOL']
    keys2 = ['baseline','stoch1_long','stoch2_long','stoch3_long',\
            'stoch4_long','stoch5_long','SPOL']   
    keys3 = ['baseline','theta_pert1','theta_pert2','theta_pert3',\
            'theta_pert4','theta_pert5','SPOL']      

    
    colors1 = ['black',\
               'lightsalmon',\
               'salmon',\
               'red',\
               'darkred',\
               'firebrick','grey']
    
    colors2 = ['black',\
         'powderblue',\
         'deepskyblue',\
         'dodgerblue',\
         'blue',\
         'navy','grey']
    
    colors3 = ['black',\
             'orchid',\
             'fuchsia',\
             'mediumorchid',\
             'purple',\
             'darkorchid','grey']
    

    dum = 0
    #for key,val in area_compdbz.items():
    for key in keys1:
        val = echo_top_20_max[key]
        val = np.array(val)
        if key != 'SPOL':
            val = val/1.e3
        # filter
        val = val[val > 6.]
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax1.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
       # tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])
        #ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])#
        dum = dum + 1    
        
    dum = 0
    #for key,val in area_compdbz.items():
    for key in keys2:
        val = echo_top_20_max[key]
        val = np.array(val)     
        if key != 'SPOL':
            val = val/1.e3
        # filter
        val = val[val > 6.]            
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax2.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax5.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
       # tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])
        #ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])#
        dum = dum + 1    
        
    dum = 0
    #for key,val in area_compdbz.items():
    for key in keys3:
        val = echo_top_20_max[key]
        val = np.array(val)
        if key != 'SPOL':
            val = val/1.e3
        # filter
        val = val[val > 6.]            
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax6.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
       # tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])
        #ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])#
        dum = dum + 1            
   

    #ax3.legend(fontsize=15,bbox_to_anchor=(2.,-0.18),ncol=4)
    ax4.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)
    ax5.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)
    ax6.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)

    ax4.set_yscale('log')
    ax5.set_yscale('log')
    ax6.set_yscale('log')
    #ax2.set_yscale('log')
    #ax4.set_yscale('log')
    #ax6.set_yscale('log')
    #ax8.set_yscale('log')
    #ax1.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax2.set_title('PDF, $Log-y$',fontsize=25)
    #ax3.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax4.set_title('Histogram, $Log-y$',fontsize=25)
    #ax5.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax6.set_title('PDF, $Log-y$',fontsize=25)
    #ax7.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax8.set_title('Histogram, $Log-y$',fontsize=25)
    
    
    plt.subplots_adjust(wspace=0.3,hspace=0.35)
    
    
    plt.show()
    print(aaaaaa)
    
    
#----------------------------------
# PDF of cell max 40 dBZ ETH
#----------------------------------

cellmax_40dbz_run = 0
if cellmax_40dbz_run == 1:

    fig = plt.figure(figsize=(18,10))

    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
    for ax in axlist:
        ax.grid(which='both',ls='solid',color='grey',lw=0.5)
        ax.set_xlabel('Core Max. 40 dBZ ETH [km]',fontsize=15)
        ax.tick_params(labelsize=15)
        #ax.set_xscale('log')
        #ax.set_xlim(20,0.5e3)
        ax.set_ylabel('Frequency',fontsize=15)

    
    #bins = np.arange(0,21,1)
    bins = np.arange(0,21,1)
    midbins = []
    for ii in range(len(bins)-1):
        midbins.append(0.5*(bins[ii]+bins[ii+1]))
    midbins = np.array(midbins)

    keys1 = ['baseline','stoch1_short','stoch2_short','stoch3_short',\
            'stoch4_short','stoch5_short','SPOL']
    keys2 = ['baseline','stoch1_long','stoch2_long','stoch3_long',\
            'stoch4_long','stoch5_long','SPOL']   
    keys3 = ['baseline','theta_pert1','theta_pert2','theta_pert3',\
            'theta_pert4','theta_pert5','SPOL']      

    
    colors1 = ['black',\
               'lightsalmon',\
               'salmon',\
               'red',\
               'darkred',\
               'firebrick','grey']
    
    colors2 = ['black',\
         'powderblue',\
         'deepskyblue',\
         'dodgerblue',\
         'blue',\
         'navy','grey']
    
    colors3 = ['black',\
             'orchid',\
             'fuchsia',\
             'mediumorchid',\
             'purple',\
             'darkorchid','grey']
    

    dum = 0
    #for key,val in area_compdbz.items():
    for key in keys1:
        val = echo_top_40_max[key]
        val = np.array(val)
        if key != 'SPOL':
            val = val/1.e3
        # filter
        val = val[val > 4.]
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax1.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
       # tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])
        #ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])#
        dum = dum + 1    
        
    dum = 0
    #for key,val in area_compdbz.items():
    for key in keys2:
        val = echo_top_40_max[key]
        val = np.array(val)     
        if key != 'SPOL':
            val = val/1.e3
        # filter
        val = val[val > 4.]            
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax2.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax5.plot(midbins,tmp_hist,label=tmplabel,color=colors2[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
       # tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])
        #ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])#
        dum = dum + 1    
        
    dum = 0
    #for key,val in area_compdbz.items():
    for key in keys3:
        val = echo_top_40_max[key]
        val = np.array(val)
        if key != 'SPOL':
            val = val/1.e3
        # filter
        val = val[val > 4.]            
        # PDF
        tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=True)
        tmp_hist[tmp_hist == 0] = np.nan
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])
        ax6.plot(midbins,tmp_hist,label=tmplabel,color=colors3[dum],\
                 lw=2,ls='solid')#,marker=markers[dum])

        # Histogram
        #tmp_hist,tmp_bins =np.histogram(val,bins=bins,normed=False)
       # tmplabel=key+' N={}'.format(str(int(len(val))))
        #ax3.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])
        #ax4.plot(midbins,tmp_hist,label=tmplabel,color=colors1[dum],\
        #         lw=lws1[dum],ls=lss1[dum])#,marker=markers[dum])#
        dum = dum + 1            
   

    #ax3.legend(fontsize=15,bbox_to_anchor=(2.,-0.18),ncol=4)
    ax4.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)
    ax5.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)
    ax6.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)

    ax4.set_yscale('log')
    ax5.set_yscale('log')
    ax6.set_yscale('log')
    #ax2.set_yscale('log')
    #ax4.set_yscale('log')
    #ax6.set_yscale('log')
    #ax8.set_yscale('log')
    #ax1.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax2.set_title('PDF, $Log-y$',fontsize=25)
    #ax3.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax4.set_title('Histogram, $Log-y$',fontsize=25)
    #ax5.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax6.set_title('PDF, $Log-y$',fontsize=25)
    #ax7.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax8.set_title('Histogram, $Log-y$',fontsize=25)
    
    
    plt.subplots_adjust(wspace=0.3,hspace=0.35)
    
    
    plt.show()
    print(aaaaaa)    
    
    
    
    
    
    
area_boxplot_run = 0
if area_boxplot_run == 1:
    

    fig = plt.figure(figsize=(15,18))

    ax = fig.add_subplot(111)
    #ax2 = fig.add_subplot(111)
    #ax3 = fig.add_subplot(111)

    ax.grid(which='both',ls='solid',color='grey',lw=0.5)
    ax.set_xlabel('Cell Area [km$^{2}$]',fontsize=25)
    ax.tick_params(labelsize=25)
    ax.set_xscale('log')
    ax.set_xlim(6,1.e3)
    
    c='grey'
    
    labels = sim_names.copy()
    labels.append('SPOL')
    
    bplot=ax.boxplot([area_compdbz['baseline'],\
                    area_compdbz['stoch1_short'],\
                    area_compdbz['stoch2_short'],\
                    area_compdbz['stoch3_short'],\
                    area_compdbz['stoch4_short'],\
                    area_compdbz['stoch5_short'],\
                    area_compdbz['stoch1_long'],\
                    area_compdbz['stoch2_long'],\
                    area_compdbz['stoch3_long'],\
                    area_compdbz['stoch4_long'],\
                    area_compdbz['stoch5_long'],\
                    area_compdbz['theta_pert1'],\
                    area_compdbz['theta_pert2'],\
                    area_compdbz['theta_pert3'],\
                    area_compdbz['theta_pert4'],\
                    area_compdbz['theta_pert5'],\
                    area_compdbz['SPOL'],\
                   ],\
                  notch=True, patch_artist=True,boxprops=dict(facecolor=c, color=c),\
              vert=False)
        #ax.boxplot(les_area_height_dict[hlevs[dum]])        
        
    ax.set_yticklabels(labels,fontsize=15,rotation=45)
    ax.grid(axis='x')        
    
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)  
        
    positions = np.arange(17) + 1
    means = [np.mean(val) for kev,val in area_compdbz.items()]
    medians = [np.median(val) for kev,val in area_compdbz.items()]
    ax.plot(means,positions,'s',markersize=20,c='magenta',zorder=100,markeredgecolor='black')  
    ax.plot(medians,positions,'|',markersize=30,\
            markeredgewidth=6,c='gold',zorder=100)  
    

    ax.axvline(np.median(area_compdbz['baseline']),ls='solid',lw=10,color='black',label='BASELINE Median')
    
    ax.legend(fontsize=25,bbox_to_anchor=(0.5,-0.15),loc='lower center',ncol=3)
    
    tmpnum = []
    for key,val in area_compdbz.items():
        tmpnum.append(len(val))
        
    dum=0
    for position in positions:
        ax.text(7,position,'N='+str(tmpnum[dum]),fontsize=20,\
                va='center',fontweight='bold')
        dum+=1

    print(aaaa)

#import gc  
#gc.collect()
    
    
#----------------------------------
# CDF and histogram of cell area
#----------------------------------



area_cdf_run = 0
if area_cdf_run == 1:
    

    fig = plt.figure(figsize=(18,12))

    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
    for ax in axlist:
        ax.grid(which='both',ls='solid',color='grey',lw=0.5)
        ax.set_xlabel('Core Area [km$^{2}$]',fontsize=15)
        ax.tick_params(labelsize=15)
        ax.set_xscale('log')


    keys1 = ['baseline','stoch1_short','stoch2_short','stoch3_short',\
            'stoch4_short','stoch5_short','SPOL']
    keys2 = ['baseline','stoch1_long','stoch2_long','stoch3_long',\
            'stoch4_long','stoch5_long','SPOL']   
    keys3 = ['baseline','theta_pert1','theta_pert2','theta_pert3',\
            'theta_pert4','theta_pert5','SPOL']      

    
    colors1 = ['black',\
               'lightsalmon',\
               'salmon',\
               'red',\
               'darkred',\
               'firebrick','grey']
    
    colors2 = ['black',\
         'powderblue',\
         'deepskyblue',\
         'dodgerblue',\
         'blue',\
         'navy','grey']
    
    colors3 = ['black',\
             'orchid',\
             'fuchsia',\
             'mediumorchid',\
             'purple',\
             'darkorchid','grey']
    

    dum = 0
    for key in keys1:
        val = area[key]
        val = np.array(val)
        val_10dbz = np.array(max_echo_top_10[key])
        if key == 'SPOL':
            val_10dbz = val_10dbz*1000.
        #tmpid = np.where(val_10dbz >= 8000.)
        #tmpid = np.where((val_10dbz >= 6000.) & (val_10dbz < 8000.))
        
        #val = val[tmpid]
        #val = val[val > 9.]
        
        # CDF
        X2 = np.sort(val)
        N = len(X2)
        F2 = np.array(range(N))/float(N)
        
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax1.plot(X2,F2,label=tmplabel,color=colors1[dum],lw=2,ls='solid')
        # CDF contribution to total
        ax4.plot(X2,np.cumsum(X2)/np.sum(X2),label=tmplabel,\
                 color=colors1[dum],lw=2,ls='solid')
        dum = dum + 1    
        
        
    dum = 0
    for key in keys2:
        val = area[key]
        val = np.array(val) 
        val_10dbz = np.array(max_echo_top_10[key])
        if key == 'SPOL':
            val_10dbz = val_10dbz*1000.       
        #tmpid = np.where(val_10dbz >= 8000.)
        #tmpid = np.where((val_10dbz >= 6000.) & (val_10dbz < 8000.))
        #val = val[tmpid]        
        #val = val[val > 9.]
        
        
        # CDF
        X2 = np.sort(val)
        N = len(X2)
        F2 = np.array(range(N))/float(N)
        
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax2.plot(X2,F2,label=tmplabel,color=colors2[dum],lw=2,ls='solid')
        # CDF contribution to total
        ax5.plot(X2,np.cumsum(X2)/np.sum(X2),label=tmplabel,\
                 color=colors2[dum],lw=2,ls='solid')
        dum = dum + 1   
        
        
        
    dum = 0
    for key in keys3:
        val = area[key]
        val = np.array(val)
        val_10dbz = np.array(max_echo_top_10[key])
        #print(np.max(val_10dbz))
        if key == 'SPOL':
            val_10dbz = val_10dbz*1000.
        #tmpid = np.where(val_10dbz >= 8000.)
        #tmpid = np.where((val_10dbz >= 6000.) & (val_10dbz < 8000.))
        #val = val[tmpid]              
        #val = val[val > 9.]
        
        
        # CDF
        X2 = np.sort(val)
        N = len(X2)
        F2 = np.array(range(N))/float(N)
        
        tmplabel=key+' N={}'.format(str(int(len(val))))
        ax3.plot(X2,F2,label=tmplabel,color=colors3[dum],lw=2,ls='solid')
        # CDF contribution to total
        ax6.plot(X2,np.cumsum(X2)/np.sum(X2),label=tmplabel,\
                 color=colors3[dum],lw=2,ls='solid')
        dum = dum + 1   
        
   

    #ax3.legend(fontsize=15,bbox_to_anchor=(2.,-0.18),ncol=4)
    ax4.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)
    ax5.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)
    ax6.legend(fontsize=15,bbox_to_anchor=(0.9,-0.2),ncol=1)

    #ax4.set_yscale('log')
    #ax5.set_yscale('log')
    #ax6.set_yscale('log')
    #ax2.set_yscale('log')
    #ax4.set_yscale('log')
    #ax6.set_yscale('log')
    #ax8.set_yscale('log')
    #ax1.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax2.set_title('PDF, $Log-y$',fontsize=25)
    #ax3.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax4.set_title('Histogram, $Log-y$',fontsize=25)
    #ax5.set_title('PDF, $Lineary-y$',fontsize=25)
    #ax6.set_title('PDF, $Log-y$',fontsize=25)
    #ax7.set_title('Histogram, $Linear-y$',fontsize=25)
    #ax8.set_title('Histogram, $Log-y$',fontsize=25)
    
    #plt.suptitle('8 Dec.',fontsize=40)
    #plt.suptitle('All Times, 6 km <= 10 dBZ ETH < 8 km',fontsize=40)
    #plt.suptitle('Before 21 UTC on 7th',fontsize=40)
    
    
    
    ax1.set_ylabel('Cumulative Frequency',fontsize=15)
    ax2.set_ylabel('Cumulative Frequency',fontsize=15)
    ax3.set_ylabel('Cumulative Frequency',fontsize=15)
    ax4.set_ylabel('Contribution of Core Area to\n Total Core Area',fontsize=15)
    ax5.set_ylabel('Contribution of Core Area to\n Total Core Area',fontsize=15)    
    ax6.set_ylabel('Contribution of Core Area to\n Total Core Area',fontsize=15)    
    
    plt.subplots_adjust(wspace=0.3,hspace=0.35)
    
    
    plt.show()
    print(aaaaaa)

    
    
