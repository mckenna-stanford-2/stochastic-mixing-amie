#==========================================================================
# Title: plot_vol_acc_precip_fig_5.py
# Author: McKenna W. Stanford
# Utility: Plot time series of domain accumulated volumetric precipitation
# within SPOL Domain. Note that for comparison with SPOL, this uses
# instantaneous rain rates.
# Used to make Fig. 5 of mansuscript.
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

path = '/glade/scratch/mckenna/AMIE/amie_post/inst_rain_rates/'

# Make SPOL range ring
x0 = np.arange(-148.5,151.5,3)
y0 = np.arange(-148.5,151.5,3)
X0,Y0 = np.meshgrid(x0,y0)
radial_dist = np.sqrt(X0**2 + Y0**2)

#=========================================================
# Make dictionary to hold all time series
vol_acc_dict = {}
#=========================================================


#=========================================================
# Calculate domain-accumulated volumetric
# precipitation for 1km CG simulation
# Optional in case you've already
# calculated it.
#=========================================================


icalc_1kmcg = 1.
if icalc_1kmcg == 1.:
    print('Processing 1km_CG simulation...')

    # Load file
    pkl_file = open(path+'baseline_1km_CG_inst_rain_rates.p','rb')
    tmpfile = pickle.load(pkl_file)
    pkl_file.close()  
    lon = tmpfile['lon']
    lat = tmpfile['lat']
    time = tmpfile['time']

    if False:
        for key,val in tmpfile.items():
            print(key,np.shape(val),np.max(val),np.min(val))
        
    
    rain_rate = tmpfile['rain_rate']
    # Limit to only the center of the domain where SPOL is
    rain_rate = rain_rate[:,100-50:100+50,100-50:100+50]
    lon = lon[100-50:100+50,100-50:100+50]
    lat = lat[100-50:100+50,100-50:100+50]
    
    
    # set points outside of the SPOL 150-km radius to NaN
    tmpid = np.where(radial_dist > 150.)
    rain_rate[:,tmpid[0],tmpid[1]] = np.nan

        
    iplot = False
    if iplot:
        rain_rate_levs = 10.**np.arange(-2,2.1,0.1)
        fig = plt.figure(figsize=(10,4))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        dum_ticks = [1.e-2,1.e-1,1.e0,1.e1,1.e2]
        plot1=ax1.contourf(tmpfile['lon'],tmpfile['lat'],tmpfile['rain_rate'][-1,:,:],levels=rain_rate_levs,\
                           cmap='cividis',norm=matplotlib.colors.LogNorm(),extend='both')
        cbar1 = fig.colorbar(plot1,ax=ax1,ticks=dum_ticks)
        cbar1.ax.set_ylabel('Rain Rate [mm hr$^{-1}$]')

        plot2=ax2.contourf(lon,lat,rain_rate[-1,:,:],levels=rain_rate_levs,\
                           cmap='cividis',norm=matplotlib.colors.LogNorm(),extend='both')
        cbar2 = fig.colorbar(plot2,ax=ax2,ticks=dum_ticks)
        cbar2.ax.set_ylabel('Rain Rate [mm hr$^{-1}$]')

        ax1.set_title('Full Domain')
        ax2.set_title('SPOL-limited Domain')
        plt.show()
        plt.close()        
        
        #plt.contourf(lon,lat,rr_all[-1,:,:],levels=levs,cmap='inferno')
        #plt.contourf(np.transpose(rr[96,:,:]),levels=[0,0.01,0.05,0.1,0.5,1,5])

        
    # Sum along lon and lat axes and vonvert to volumetric precipitatio (mm km^2)
    rain_rate_sum = np.nansum(rain_rate,axis=(1,2))*9
    
    vol_acc_dict['baseline_1km_CG'] = rain_rate_sum
    vol_acc_dict['time'] = time
    
    iplot = False
    if iplot:
        Fontsize=14
        fig = plt.figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        ax.grid(ls='dotted',c='dimgrey',lw=1)
        ax.set_xlabel('Time [UTC]',fontsize=Fontsize)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Domain-Accumulated\nVolumetric Rainfall\n[x10$^{5}$ mm km$^{2}$]',fontsize=Fontsize)#,labelpad=30)
        ax.tick_params(labelsize=Fontsize)
        ax.plot(time,rain_rate_sum*1.e-5,lw=2,c='navy')
        plt.show()
        plt.close()        

        
#=========================================================
# Calculate domain-accumulated volumetric
# precipitation for SPOL.
# Optional in case you've already
# calculated it.
#=========================================================

icalc_spol = 1.
if icalc_spol == 1.:
    print('Processing SPOL...')

    spol_path = '/glade/scratch/mckenna/AMIE/amie_post/spol/'
    # Load file
    pkl_file = open(spol_path+'spol_coarse_grained.p','rb')
    tmpfile = pickle.load(pkl_file)
    pkl_file.close()
    
    # Grab 3 km coarse-grained dictionary
    tmpfile = tmpfile['3km']
    
    if False:
        for key,val in tmpfile.items():
            print(key,np.shape(val),np.max(val),np.min(val))
    
    
    lon = tmpfile['lon']
    lat = tmpfile['lat']
    time = tmpfile['time']
    rain_rate = tmpfile['rain_rate']
    rain_rate_min = tmpfile['rain_rate_min']
    rain_rate_max = tmpfile['rain_rate_max']
    
       
    iplot = False
    if iplot:
        rain_rate_levs = 10.**np.arange(-2,2.1,0.1)
        fig = plt.figure(figsize=(6,4))
        ax1 = fig.add_subplot(111)
        dum_ticks = [1.e-2,1.e-1,1.e0,1.e1,1.e2]
        plot1=ax1.contourf(lon,lat,rain_rate[-1,:,:],levels=rain_rate_levs,\
                           cmap='cividis',norm=matplotlib.colors.LogNorm(),extend='both')
        cbar1 = fig.colorbar(plot1,ax=ax1,ticks=dum_ticks)
        cbar1.ax.set_ylabel('Rain Rate [mm hr$^{-1}$]')

        plt.show()
        plt.close()        
        
    # Sum along lon and lat axes and vonvert to volumetric precipitatio (mm km^2)
    rain_rate_sum = np.nansum(rain_rate,axis=(1,2))*9
    rain_rate_max_sum = np.nansum(rain_rate_max,axis=(1,2))*9
    rain_rate_min_sum = np.nansum(rain_rate_min,axis=(1,2))*9
    
    vol_acc_dict['spol'] = rain_rate_sum
    vol_acc_dict['spol_max'] = rain_rate_max_sum
    vol_acc_dict['spol_min'] = rain_rate_min_sum
    
    iplot = False
    if iplot:
        Fontsize=14
        fig = plt.figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        ax.grid(ls='dotted',c='dimgrey',lw=1)
        ax.set_xlabel('Time [UTC]',fontsize=Fontsize)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_ylabel('Domain-Accumulated\nVolumetric Rainfall\n[x10$^{5}$ mm km$^{2}$]',fontsize=Fontsize)
        ax.tick_params(labelsize=Fontsize)
        ax.plot(time,rain_rate_sum*1.e-5,lw=2,c='navy',label='Best')
        ax.plot(time,rain_rate_min_sum*1.e-5,lw=2,c='deepskyblue',ls='dashed',label='Min')
        ax.plot(time,rain_rate_max_sum*1.e-5,lw=2,c='deepskyblue',ls='dashed',label='Max')
        ax.legend(loc='upper left',fontsize=Fontsize)
        plt.show()
        plt.close()  
        
        
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
        pkl_file = open(path+sim+'_inst_rain_rates.p','rb')
        tmpfile = pickle.load(pkl_file)
        pkl_file.close()  
        lon = tmpfile['lon']
        lat = tmpfile['lat']
        time = tmpfile['time']

        if False:
            for key,val in tmpfile.items():
                print(key,np.shape(val),np.max(val),np.min(val))


        rain_rate = tmpfile['rain_rate']
        # Limit to only the center of the domain where SPOL is
        rain_rate = rain_rate[:,200-50:200+50,200-50:200+50]
        lon = lon[200-50:200+50,200-50:200+50]
        lat = lat[200-50:200+50,200-50:200+50]


        # set points outside of the SPOL 150-km radius to NaN
        tmpid = np.where(radial_dist > 150.)
        rain_rate[:,tmpid[0],tmpid[1]] = np.nan


        iplot = False
        if iplot:
            rain_rate_levs = 10.**np.arange(-2,2.1,0.1)
            fig = plt.figure(figsize=(10,4))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            dum_ticks = [1.e-2,1.e-1,1.e0,1.e1,1.e2]
            plot1=ax1.contourf(tmpfile['lon'],tmpfile['lat'],tmpfile['rain_rate'][-1,:,:],levels=rain_rate_levs,\
                               cmap='cividis',norm=matplotlib.colors.LogNorm(),extend='both')
            cbar1 = fig.colorbar(plot1,ax=ax1,ticks=dum_ticks)
            cbar1.ax.set_ylabel('Rain Rate [mm hr$^{-1}$]')

            plot2=ax2.contourf(lon,lat,rain_rate[-1,:,:],levels=rain_rate_levs,\
                               cmap='cividis',norm=matplotlib.colors.LogNorm(),extend='both')
            cbar2 = fig.colorbar(plot2,ax=ax2,ticks=dum_ticks)
            cbar2.ax.set_ylabel('Rain Rate [mm hr$^{-1}$]')

            ax1.set_title('Full Domain')
            ax2.set_title('SPOL-limited Domain')
            plt.show()
            plt.close()        

        # Sum along lon and lat axes and vonvert to volumetric precipitatio (mm km^2)
        rain_rate_sum = np.nansum(rain_rate,axis=(1,2))*9
        vol_acc_dict[sim] = rain_rate_sum

        iplot = False
        if iplot:
            Fontsize=14
            fig = plt.figure(figsize=(10,4))
            ax = fig.add_subplot(111)
            ax.grid(ls='dotted',c='dimgrey',lw=1)
            ax.set_xlabel('Time [UTC]',fontsize=Fontsize)
            ax.xaxis.set_major_formatter(dfmt)
            ax.set_ylabel('Domain-Accumulated\nVolumetric Rainfall\n[x10$^{5}$ mm km$^{2}$]',fontsize=Fontsize)#,labelpad=30)
            ax.tick_params(labelsize=Fontsize)
            ax.plot(time,rain_rate_sum*1.e-5,lw=2,c='navy')
            plt.show()
            plt.close()
            
        
print('Diagnostics')   
for key,val in vol_acc_dict.items():
    print(key,np.shape(val),np.max(val),np.min(val))
        
        
#=========================================================       
# Perform running means
#=========================================================
run_mean_dict = {}
run_mean_dict['time'] = vol_acc_dict['time']
irun_mean = 1.
if irun_mean == 1.:
    N=4
    for key,val in vol_acc_dict.items():
        if key == 'time':
            continue
        tmpval = val.copy()
        tmpval_rm = running_mean(tmpval,N)
        run_mean_dict[key] = tmpval_rm

#=============================================================
# Plot Fig. 5 for manuscript
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
          'baseline_1km_CG',\
          'SPOL']
   

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
        'no_mixing','4x','baseline_1km_CG']
lss3 = ['solid','solid','solid','solid','solid','solid','dotted','dashed','-.']

mult_fac = 1.e-5


#--------------------------------
# Start figure
#--------------------------------

fig = plt.figure(figsize=(8.3,16))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
Fontsize=20
fac = 1.25

axlist = [ax1,ax2,ax3]
for ax in axlist:
    ax.grid(ls='dotted',lw=1,c='dimgrey')
    ax.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize*fac)
    ax.set_xlim(run_mean_dict['time'].min(),run_mean_dict['time'].max())
    ax.xaxis.set_major_formatter(dfmt)
    ax.set_ylabel('Domain-Accumulated\nVolumetric Rainfall\n[x10$^{5}$ mm km$^{2}$]',fontsize=Fontsize*fac)
    ax.tick_params(labelsize=Fontsize*fac)
    ax.xaxis.set_ticks(time[::24])
    ax.plot(run_mean_dict['time'],run_mean_dict['spol']*mult_fac,lw=3,c='grey',label='SPOL')
    ax.fill_between(run_mean_dict['time'],\
                    run_mean_dict['spol']*mult_fac,\
                    run_mean_dict['spol_max']*mult_fac,color='grey',alpha=0.25)
    ax.fill_between(run_mean_dict['time'],\
                    run_mean_dict['spol_min']*mult_fac,\
                    run_mean_dict['spol']*mult_fac,color='grey',alpha=0.25)    
    
    # STOCH_SHORT Simulations
    dum = 0
    for sim in sims1:
        ax1.plot(run_mean_dict['time'],run_mean_dict[sim]*mult_fac,color=colors1[dum],lw=lws1[dum],ls=lss1[dum])
        dum+=1
        
    # STOCH_LONG Simulations
    dum = 0
    for sim in sims2:
        ax2.plot(run_mean_dict['time'],run_mean_dict[sim]*mult_fac,color=colors2[dum],lw=lws2[dum],ls=lss2[dum])
        dum+=1
        
    # Perturbed parameter and theta_pert simulations
    dum = 0
    for sim in sims3:
        ax3.plot(run_mean_dict['time'],run_mean_dict[sim]*mult_fac,color=colors3[dum],lw=lws3[dum],ls=lss3[dum])
        dum+=1
    
plt.text(0.5,1.07,'STOCH\_SHORT',fontsize=Fontsize*2,fontweight='bold',\
         bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
         rotation=0,ha='center',va='bottom')
plt.text(0.5,1.07,'STOCH\_LONG',fontsize=Fontsize*2,fontweight='bold',\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
         rotation=0,ha='center',va='bottom')
plt.text(0.5,1.07,'$\\theta$-pert',fontsize=Fontsize*2,fontweight='bold',\
         bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
         rotation=0,ha='center',va='bottom')  

ax1.text(0.02,0.77,'(a)',fontsize=Fontsize*2.5,fontweight='bold',transform=ax1.transAxes)
ax2.text(0.02,0.77,'(b)',fontsize=Fontsize*2.5,fontweight='bold',transform=ax2.transAxes)
ax3.text(0.02,0.77,'(c)',fontsize=Fontsize*2.5,fontweight='bold',transform=ax3.transAxes)


custom_lines_1 = [Line2D([0],[0],color = 'black',lw=3),\
                Line2D([0],[0],color='grey',lw=3),\
                Line2D([0],[0],color='darkorange',ls='-.',lw=3),\
               ]

#ax1.legend(custom_lines_1,['BASELINE','SPOL',],fontsize=Fontsize,\
#                loc='upper left',ncol=1,framealpha=1)
#ax2.legend(custom_lines_1,['BASELINE','SPOL',],fontsize=Fontsize,\
#                loc='upper left',ncol=1,framealpha=1)   

custom_lines_2 = [Line2D([0],[0],color = 'black',lw=3),\
                Line2D([0],[0],color='grey',lw=3),\
                Line2D([0],[0],color='black',ls='dotted',lw=3),\
                Line2D([0],[0],color='black',ls='dashed',lw=3),\
                Line2D([0],[0],color='darkorange',ls='-.',lw=3),\
               ]
ax3.legend(custom_lines_2,['BASELINE','SPOL','NO\_MIXING','4X','BASELINE\_1KM\_CG'],fontsize=Fontsize,\
                loc='lower center',ncol=2,framealpha=1,bbox_to_anchor=(0.5,-0.85)) 



plt.subplots_adjust(hspace=0.6)

outfile = 'fig_05.png'
plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
#plt.show()
plt.close()
    

    
    
    
    
    
    
    
    
    
    
    
    
#===========================================
# Domain-accumulated volumetric precip
# at end of simulation and relative differences
# (for table)
#===========================================

# limit domain-total calculations to only after 12Z

tmp_spol_vol_acc = spol_vol_acc[48:]
tmp_spol_min_vol_acc = spol_min_vol_acc[48:]
tmp_spol_max_vol_acc = spol_max_vol_acc[48:]

tmp_spol_vol_acc_tot = np.cumsum(spol_vol_acc)
tmp_spol_min_vol_acc_tot = np.cumsum(spol_min_vol_acc)
tmp_spol_max_vol_acc_tot = np.cumsum(spol_max_vol_acc)

tmp_rain_vol_acc_tot = {}

for key,val in rain_vol_acc.items():
    tmpval = val[48:]
    tmp_rain_vol_acc_tot[key] = np.cumsum(tmpval)

last_value = {}

for key,val in tmp_rain_vol_acc_tot.items():
    last_value[key] = val[-1]


last_value['SPOL'] = tmp_spol_vol_acc_tot[-1] 
last_value['SPOL_min'] = tmp_spol_min_vol_acc_tot[-1] 
last_value['SPOL_max'] = tmp_spol_max_vol_acc_tot[-1] 

for key,val in last_value.items():
    print(key,val)
    val = val/1.e6
    last_value[key] = val


# now compute relative difference from observations
rel_diffs = {}
for key,val in last_value.items():
    if (key != 'SPOL') & (key != 'SPOL_min') & (key != 'SPOL_max'):
        rel_diffs[key] = (val-last_value['SPOL'])/(last_value['SPOL'])*100.

for key,val in rel_diffs.items():
    print(key,val)
    
print(aaaaa)    
    
