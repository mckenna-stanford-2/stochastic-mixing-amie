#==========================================================================
# Title: plot_conv_strat_rain_rates_fig7.py
# Author: McKenna W. Stanford
# Utility:# Utility: Plots time series of convective and stratiform
# rain rates limited to the SPOL domain.
# Used as Fig. 7 in manuscript
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
    
    
path = '/glade/scratch/mckenna/AMIE/amie_post/'

# Define SPOL range ring
x0 = np.arange(-148.5,151.5,3)
y0 = np.arange(-148.5,151.5,3)
X0,Y0 = np.meshgrid(x0,y0)
radial_dist = np.sqrt(X0**2 + Y0**2)



#=========================================================
# SPOL convective and stratiform IDs
#=========================================================
pkl_file = open(path+'conv_strat_id/spol_conv_strat_id.p','rb')
spol_conv_strat_dict = pickle.load(pkl_file)
pkl_file.close()

spol_con3_mask = spol_conv_strat_dict['con3_mask']

#=========================================================
# SPOL rain rates
#=========================================================
print('Beginning SPOL calculations.')

pkl_file = open(path+'spol/spol_coarse_grained.p','rb')
tmpfile = pickle.load(pkl_file)['3km']
pkl_file.close()   
    
spol_rain_rate = tmpfile['rain_rate']
spol_rain_rate_min = tmpfile['rain_rate_min']
spol_rain_rate_max = tmpfile['rain_rate_max']
spol_lon = tmpfile['lon']
spol_lat = tmpfile['lat']
spol_time = tmpfile['time']
nt = len(spol_time)
    

spol_conv_rain_rate_mean = []
spol_conv_rain_rate_max_mean = []
spol_conv_rain_rate_min_mean = []
spol_strat_rain_rate_mean = []
spol_strat_rain_rate_max_mean = []
spol_strat_rain_rate_min_mean = []
    
# loop through times and calculation
# domain-mean convective and stratiform
# rain rates
    
for tt in range(nt):
    tmp_spol_rain_rate = spol_rain_rate[tt,:,:]
    tmp_spol_rain_rate_min = spol_rain_rate_min[tt,:,:]
    tmp_spol_rain_rate_max = spol_rain_rate_max[tt,:,:]
        
    tmp_spol_conv_strat_id = spol_con3_mask[tt,:,:]
    conv_id = np.where(tmp_spol_conv_strat_id == 1)
    strat_id = np.where(tmp_spol_conv_strat_id == 2)
    
    tmp_conv_rain_rate_mean = np.mean(tmp_spol_rain_rate[conv_id])
    tmp_conv_rain_rate_max_mean = np.mean(tmp_spol_rain_rate_max[conv_id])
    tmp_conv_rain_rate_min_mean = np.mean(tmp_spol_rain_rate_min[conv_id])
    tmp_strat_rain_rate_mean = np.mean(tmp_spol_rain_rate[strat_id])
    tmp_strat_rain_rate_max_mean = np.mean(tmp_spol_rain_rate_max[strat_id])
    tmp_strat_rain_rate_min_mean = np.mean(tmp_spol_rain_rate_min[strat_id])
    
    
    spol_conv_rain_rate_mean.append(tmp_conv_rain_rate_mean)
    spol_conv_rain_rate_max_mean.append(tmp_conv_rain_rate_max_mean)
    spol_conv_rain_rate_min_mean.append(tmp_conv_rain_rate_min_mean)
    spol_strat_rain_rate_mean.append(tmp_strat_rain_rate_mean)
    spol_strat_rain_rate_max_mean.append(tmp_strat_rain_rate_max_mean)
    spol_strat_rain_rate_min_mean.append(tmp_strat_rain_rate_min_mean)
    

spol_conv_rain_rate_mean = np.array(spol_conv_rain_rate_mean)
spol_conv_rain_rate_max_mean = np.array(spol_conv_rain_rate_max_mean)
spol_conv_rain_rate_min_mean = np.array(spol_conv_rain_rate_min_mean)
spol_strat_rain_rate_mean = np.array(spol_strat_rain_rate_mean)
spol_strat_rain_rate_max_mean = np.array(spol_strat_rain_rate_max_mean)
spol_strat_rain_rate_min_mean = np.array(spol_strat_rain_rate_min_mean)
print('Completed SPOL calculations.')

    
iplot = False
if iplot:
    Fontsize=14
    fig = plt.figure(figsize=(12,4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    axlist = [ax1,ax2]
    for ax in axlist:
        ax.tick_params(labelsize=Fontsize)
        ax.grid(ls='dotted',lw=1,c='dimgrey')
        ax.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize)
        ax.xaxis.set_major_formatter(dfmt)
    # Convective
    ax1.plot(spol_time,spol_conv_rain_rate_mean,lw=2,c='navy')
    ax1.plot(spol_time,spol_conv_rain_rate_min_mean,lw=2,c='deepskyblue',ls='dotted')
    ax1.plot(spol_time,spol_conv_rain_rate_max_mean,lw=2,c='deepskyblue',ls='dotted')
    ax1.xaxis.set_ticks(spol_time[::24])
    ax1.set_ylabel('Convective Rain Rate [mm hr$^{-1}$]',fontsize=Fontsize)
    # Stratiform
    ax2.plot(spol_time,spol_strat_rain_rate_mean,lw=2,c='navy')
    ax2.plot(spol_time,spol_strat_rain_rate_min_mean,lw=2,c='deepskyblue',ls='dotted')
    ax2.plot(spol_time,spol_strat_rain_rate_max_mean,lw=2,c='deepskyblue',ls='dotted')
    ax2.xaxis.set_ticks(spol_time[::24])
    ax2.set_ylabel('Stratiform Rain Rate [mm hr$^{-1}$]',fontsize=Fontsize)
    
    plt.show()
    plt.close()


#=========================================================
# baseline_1km convective and stratiform IDs
#=========================================================
pkl_file = open(path+'conv_strat_id/baseline_1km_CG_conv_strat_id.p','rb')
onekm_conv_strat_dict = pickle.load(pkl_file)
pkl_file.close()

onekm_con3_mask = onekm_conv_strat_dict['con3_mask'][:,100-50:100+50,100-50:100+50]
# set points outside of 150-km radius to NaN
tmpid = np.where(radial_dist > 150.)
onekm_con3_mask[:,tmpid[0],tmpid[1]] = np.nan
#=======================================
# baseline_1km rain rates
#=======================================
print('Beginning 1-km calculations.')


pkl_file = open(path+'inst_rain_rates/baseline_1km_CG_inst_rain_rates.p','rb')
tmpfile = pickle.load(pkl_file)
pkl_file.close()   
    
onekm_rain_rate = tmpfile['rain_rate'][:,100-50:100+50,100-50:100+50]
onekm_time = tmpfile['time']
nt = len(onekm_time)
    
# set points outside of 150-km radius to NaN
tmpid = np.where(radial_dist > 150.)
onekm_rain_rate[:,tmpid[0],tmpid[1]] = np.nan
    
    
onekm_conv_rain_rate_mean = []
onekm_strat_rain_rate_mean = []
    
    
for tt in range(nt):
    tmp_onekm_rain_rate = onekm_rain_rate[tt,:,:]
        
    tmp_onekm_conv_strat_id = onekm_con3_mask[tt,:,:]
    conv_id = np.where(tmp_onekm_conv_strat_id == 1)
    strat_id = np.where(tmp_onekm_conv_strat_id == 2)
    
    tmp_conv_rain_rate_mean = np.mean(tmp_onekm_rain_rate[conv_id])
    tmp_strat_rain_rate_mean = np.mean(tmp_onekm_rain_rate[strat_id])
    
    
    onekm_conv_rain_rate_mean.append(tmp_conv_rain_rate_mean)
    onekm_strat_rain_rate_mean.append(tmp_strat_rain_rate_mean)
    

onekm_conv_rain_rate_mean = np.array(onekm_conv_rain_rate_mean)
onekm_strat_rain_rate_mean = np.array(onekm_strat_rain_rate_mean)
print('Completed baseline_1km calculations.')

    
iplot = True
if iplot:
    Fontsize=14
    fig = plt.figure(figsize=(12,4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    axlist = [ax1,ax2]
    for ax in axlist:
        ax.tick_params(labelsize=Fontsize)
        ax.grid(ls='dotted',lw=1,c='dimgrey')
        ax.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize)
        ax.xaxis.set_major_formatter(dfmt)
    # Convective
    ax1.plot(onekm_time,onekm_conv_rain_rate_mean,lw=2,c='navy')
    ax1.xaxis.set_ticks(onekm_time[::24])
    ax1.set_ylabel('Convective Rain Rate [mm hr$^{-1}$]',fontsize=Fontsize)
    # Stratiform
    ax2.plot(spol_time,onekm_strat_rain_rate_mean,lw=2,c='navy')
    ax2.xaxis.set_ticks(onekm_time[::24])
    ax2.set_ylabel('Stratiform Rain Rate [mm hr$^{-1}$]',fontsize=Fontsize)
    
    plt.show()
    plt.close()



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

    con3_mask = conv_strat_dict['con3_mask'][:,200-50:200+50,200-50:200+50]
    # set points outside of 150-km radius to NaN
    tmpid = np.where(radial_dist > 150.)
    con3_mask[:,tmpid[0],tmpid[1]] = np.nan

    #=======================================
    # rain rates
    #=======================================

    pkl_file = open(path+'inst_rain_rates/{}_inst_rain_rates.p'.format(sim),'rb')
    tmpfile = pickle.load(pkl_file)
    pkl_file.close()   

    rain_rate = tmpfile['rain_rate'][:,200-50:200+50,200-50:200+50]
    time = tmpfile['time']
    nt = len(time)
    
    # set points outside of 150-km radius to NaN
    tmpid = np.where(radial_dist > 150.)
    rain_rate[:,tmpid[0],tmpid[1]] = np.nan
    

    conv_rain_rate_mean = []
    strat_rain_rate_mean = []


    for tt in range(nt):
        tmp_rain_rate = rain_rate[tt,:,:]

        tmp_conv_strat_id = con3_mask[tt,:,:]
        conv_id = np.where(tmp_conv_strat_id == 1)
        strat_id = np.where(tmp_conv_strat_id == 2)

        tmp_conv_rain_rate_mean = np.mean(tmp_rain_rate[conv_id])
        tmp_strat_rain_rate_mean = np.mean(tmp_rain_rate[strat_id])


        conv_rain_rate_mean.append(tmp_conv_rain_rate_mean)
        strat_rain_rate_mean.append(tmp_strat_rain_rate_mean)


    conv_rain_rate_mean = np.array(conv_rain_rate_mean)
    strat_rain_rate_mean = np.array(strat_rain_rate_mean)    
    
    conv_rain_rate_mean_dict[sim] = conv_rain_rate_mean
    strat_rain_rate_mean_dict[sim] = strat_rain_rate_mean
    
    
    
    
    
# Add SPOL and baseline_1km sims to dictionary   
conv_rain_rate_mean_dict['spol'] = spol_conv_rain_rate_mean 
conv_rain_rate_mean_dict['spol_max'] = spol_conv_rain_rate_max_mean 
conv_rain_rate_mean_dict['spol_min'] = spol_conv_rain_rate_min_mean 
strat_rain_rate_mean_dict['spol'] = spol_strat_rain_rate_mean 
strat_rain_rate_mean_dict['spol_max'] = spol_strat_rain_rate_max_mean 
strat_rain_rate_mean_dict['spol_min'] = spol_strat_rain_rate_min_mean 
conv_rain_rate_mean_dict['baseline_1km'] = onekm_conv_rain_rate_mean 
strat_rain_rate_mean_dict['baseline_1km'] = onekm_strat_rain_rate_mean 


#=========================================
# perform running mean
#=========================================
conv_run_mean_dict = {}
strat_run_mean_dict = {}

irun_mean = 1
if irun_mean == 1:
    N=4
    # Convective
    for key,val in conv_rain_rate_mean_dict.items():
        tmpval = val.copy()
        tmpval_rm = running_mean(tmpval,N)
        conv_run_mean_dict[key] = tmpval_rm
    # Stratiform
    for key,val in strat_rain_rate_mean_dict.items():
        tmpval = val.copy()
        tmpval_rm = running_mean(tmpval,N)
        strat_run_mean_dict[key] = tmpval_rm        
        
conv_run_mean_dict['time'] = time
strat_run_mean_dict['time'] = time
    
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


#=============================================================
# Plot Fig. 7 for manuscript
#=============================================================

titles = ['baseline',\
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
          'stoch5_long',\
          'theta_pert1',\
          'theta_pert2',\
          'theta_pert3',\
          'theta_pert4',\
          'theta_pert5',\
          'no_mixing',\
          '4x',\
          'baseline_1km',\
          'spol',\
          'spol_max',\
          'spol_min']

    
colors1 = ['black',\
     'lightsalmon',\
     'salmon',\
     'red',\
     'darkred',\
     'firebrick']

sims1 = ['baseline','stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']

lws1 = [3,1.5,1.5,1.5,1.5,1.5]
lss1 = ['solid','solid','solid','solid','solid','solid']

colors2 = ['black',\
         'powderblue',\
         'deepskyblue',\
         'dodgerblue',\
         'blue',\
         'navy']

lws2 = [3,1.5,1.5,1.5,1.5,1.5]
sims2 = ['baseline','stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
lss2 = ['solid','solid','solid','solid','solid','solid']

colors3 = ['black',\
         'orchid',\
         'fuchsia',\
         'mediumorchid',\
         'purple',\
         'darkorchid',\
         'black',\
         'black',\
         'darkorange',\
         ] 
    
lws3 = [3,1.5,1.5,1.5,1.5,1.5,3,3,3]
sims3 = ['baseline','theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5',\
        'no_mixing','4x','baseline_1km']
lss3 = ['solid','solid','solid','solid','solid','solid','dotted','dashed','-.']
   
    
#--------------------------------
# Start figure
#--------------------------------

fig = plt.figure(figsize=(12,10))
ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(322)
ax3 = fig.add_subplot(323)
ax4 = fig.add_subplot(324)
ax5 = fig.add_subplot(325)
ax6 = fig.add_subplot(326)
Fontsize=16


# Convective
axlist = [ax1,ax3,ax5]
for ax in axlist:
    ax.grid(ls='dotted',lw=1,c='dimgrey')
    ax.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize)
    ax.xaxis.set_major_formatter(dfmt)
    ax.set_ylabel('Mean Convective\nRain Rate [mm hr$^{-1}$]',fontsize=Fontsize)
    ax.tick_params(labelsize=Fontsize)
    ax.xaxis.set_ticks(time[::24])
    ax.plot(time,conv_run_mean_dict['spol'],lw=3,c='grey',label='SPOL')
    ax.fill_between(time,conv_run_mean_dict['spol'],\
        conv_run_mean_dict['spol_max'],color='grey',alpha=0.25)
    ax.fill_between(time,conv_run_mean_dict['spol_min'],\
        conv_run_mean_dict['spol'],color='grey',alpha=0.25)
    ax.set_ylim(0,20)
    
# Stratiform    
axlist = [ax2,ax4,ax6]
for ax in axlist:
    ax.grid(ls='dotted',lw=1,c='dimgrey')
    ax.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize)
    ax.xaxis.set_major_formatter(dfmt)
    ax.set_ylabel('Mean Stratiform\nRain Rate [mm hr$^{-1}$]',fontsize=Fontsize)
    ax.tick_params(labelsize=Fontsize)
    ax.xaxis.set_ticks(time[::24])
    ax.plot(time,strat_run_mean_dict['spol'],lw=3,c='grey',label='SPOL')
    ax.fill_between(time,strat_run_mean_dict['spol'],\
        strat_run_mean_dict['spol_max'],color='grey',alpha=0.25)
    ax.fill_between(time,strat_run_mean_dict['spol_min'],\
        strat_run_mean_dict['spol'],color='grey',alpha=0.25)
    ax.set_ylim(0,1.1)
    
    
# STOCH_SHORT Simulations
dum = 0
for sim in sims1:
    ax1.plot(time,conv_run_mean_dict[sim],color=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    ax2.plot(time,strat_run_mean_dict[sim],color=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    dum+=1
        
# STOCH_LONG Simulations
dum = 0
for sim in sims2:
    ax3.plot(time,conv_run_mean_dict[sim],color=colors2[dum],lw=lws2[dum],ls=lss2[dum])
    ax4.plot(time,strat_run_mean_dict[sim],color=colors2[dum],lw=lws2[dum],ls=lss2[dum])
    dum+=1        
        
# Perturbed parameter and theta_pert simulations
dum = 0
for sim in sims3:
    ax5.plot(time,conv_run_mean_dict[sim],color=colors3[dum],lw=lws3[dum],ls=lss3[dum])
    ax6.plot(time,strat_run_mean_dict[sim],color=colors3[dum],lw=lws3[dum],ls=lss3[dum])
    dum+=1      
    
dumx = -0.35
dumy = 0.5
plt.text(dumx,dumy,'STOCH\_SHORT',fontsize=Fontsize*1.5,fontweight='bold',\
         bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
         rotation=90,va='center')
plt.text(dumx,dumy,'STOCH\_LONG',fontsize=Fontsize*1.5,fontweight='bold',\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax3.transAxes,\
         rotation=90,va='center')
plt.text(dumx,dumy,'$\\theta$-pert',fontsize=Fontsize*1.5,fontweight='bold',\
         bbox=dict(facecolor='purple', alpha=0.25),transform=ax5.transAxes,\
         rotation=90,va='center')  

labs = ['(a)','(b)','(c)','(d)','(e)','(f)']
axlist4 = [ax1,ax2,ax3,ax4,ax5,ax6]
dum = 0
for ax in axlist4:
    ax.text(-0.2,1.125,labs[dum],transform=ax.transAxes,fontweight='bold',fontsize=Fontsize*2.5)
    ax.set_xlim(time.min(),time.max())
    dum+=1
    

custom_lines = [Line2D([0],[0],color = 'black',lw=3),\
                Line2D([0],[0],color='grey',lw=3),\
                Line2D([0],[0],color='black',ls='dotted',lw=3),\
                Line2D([0],[0],color='black',ls='dashed',lw=3),\
                Line2D([0],[0],color='darkorange',ls='-.',lw=3),\
               ]

ax5.legend(custom_lines,['BASELINE','SPOL','NO\_MIXING','4X','BASELINE\_1KM\_CG'],\
           fontsize=Fontsize,\
            loc='lower center',bbox_to_anchor=(1.2,-0.8),ncol=3,framealpha=0)

ax1.text(0.5,1.03,'Convective Rain Rate',fontsize=Fontsize*1.5,\
         va='bottom',ha='center',transform=ax1.transAxes)
ax2.text(0.5,1.03,'Stratiform Rain Rate',fontsize=Fontsize*1.5,\
         va='bottom',ha='center',transform=ax2.transAxes)

plt.subplots_adjust(wspace=0.3,hspace=0.4)

outfile = 'fig_07.png'
plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
#plt.show()
plt.close()
    



print(aaaaaa)

