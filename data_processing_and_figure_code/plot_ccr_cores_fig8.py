#==========================================================================
# Title: plot_ccr_cores_fig8.py
# Author: McKenna W. Stanford
# Utility: Computes contiguous convective reflectivity cores and their
# properties only within the SPOL domain.
# Makes Fig. 8 of the manuscript.
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
import cartopy.crs as ccrs
from matplotlib.cm import get_cmap
import matplotlib.dates as mdates
from matplotlib.lines import Line2D

#------------------------------------------
# Functions
#------------------------------------------
# make function that finds the nearest
# element in array to the given value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx
#running mean
def running_mean(x,N):
    x_padded = np.pad(x, (N//2, N-1-N//2), mode='edge')
    x_smooth = np.convolve(x_padded, np.ones((N,))/N, mode='valid')
    return x_smooth
#------------------------------------------
# Parameters
#------------------------------------------
plt.rc('text',usetex=True)
dfmt = mdates.DateFormatter('%d-%H')

# Define SPOL range ring
wrfx0 = np.arange(-148.5,151.5,3)
wrfy0 = np.arange(-148.5,151.5,3)
tmpmeshx,tmpmeshy = np.meshgrid(wrfx0,wrfy0)
#calculate radial distance from origin
radial_dist = np.sqrt(tmpmeshx**2. + tmpmeshy**2.)

sims = ['spol',\
        'baseline_1km_CG',\
             'baseline',\
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


ccr_core_dict = {}

num_sims = len(sims)

path = '/glade/scratch/mckenna/AMIE/amie_post/conv_strat_id/'

for sim_ii in range(num_sims):
    print('Simulation:',sims[sim_ii])
    ccr_core_dict[sims[sim_ii]] = {}
    
    pkl_file = open(path+sims[sim_ii]+'_conv_strat_id.p','rb')
    tmpdict = pickle.load(pkl_file)
    
    
    con3_mask = tmpdict['con3_mask']
    lon = tmpdict['lon']
    lat = tmpdict['lat']
    time = tmpdict['time']
    nt = len(time)
    
    
    # Change stratiform indices (2) to 0
    con3_mask[con3_mask == 2.] = 0.
    
    # Limit to SPOL domain
    if sims[sim_ii] != 'spol':
        if sims[sim_ii] == 'baseline_1km_CG':
            midpoint = 100
        else:
            midpoint = 200
        halfrange = 50

        lon = lon[midpoint-halfrange:midpoint+halfrange,midpoint-halfrange:midpoint+halfrange]
        lat = lat[midpoint-halfrange:midpoint+halfrange,midpoint-halfrange:midpoint+halfrange]
        con3_mask = con3_mask[:,midpoint-halfrange:midpoint+halfrange,midpoint-halfrange:midpoint+halfrange] 
    

        # set points outside of the SPOL 150-km radius to NaN
        tmpid = np.where(radial_dist > 150.)
        con3_mask[:,tmpid[0],tmpid[1]] = np.nan
        
        
    
    #for tt in range(-3,-2):    
    for tt in range(nt):
        tmp_con3_mask = con3_mask[tt,:,:]
        regions = label(tmp_con3_mask)
        label_num = regions[1]
        regions = regions[0]
                
        # Deal with objects lying along edge of SPOL domain
        # 
        # apply a filter that says if the maximum
        # radial distance identified within the object is > 145,
        # then discard it
        for rr in range(1,np.max(regions)+1):
            rid = np.where(regions == rr)
            if np.max(radial_dist[rid]) > 145.:
                tmp_con3_mask[rid] = 0
                
        # Now relabel regions
        regions = label(tmp_con3_mask)
        label_num = regions[1]
        regions = regions[0]        
        
        iplot = False
        if iplot:
            
            fig = plt.figure(figsize=(12,5))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            Fontsize=14
            ax1.tick_params(labelsize=Fontsize)
            ax2.tick_params(labelsize=Fontsize)
            ax1.grid()
            ax2.grid()
            ax1.contourf(lon,lat,tmp_con3_mask,colors='blue',levels=[0.5,1.5],label='conv.')
            ax1.contour(lon,lat,radial_dist,levels=[149,151],colors='black',linewidths=5)   
            ax2.contour(lon,lat,radial_dist,levels=[149,151],colors='black',linewidths=5) 
            ax2.contourf(lon,lat,regions,levels=np.arange(1,label_num+1,1),cmap='nipy_spectral')
            plt.show()
            plt.close()
            
                        
        #------------------------------------
        # Now loop through cores and calculate
        # the area.
        #------------------------------------
        ccr_area = []
        for rr in range(1,label_num):
            tmpid = np.where(regions == rr)
            tmp_ccr_area = np.size(tmpid[0])*9.
            ccr_area.append(tmp_ccr_area)
        
        ccr_area = np.array(ccr_area)
        ccr_id = np.arange(1,len(ccr_area)+1,1)
        ccr_core_dict[sims[sim_ii]][str(tt+1)] = {'ccr_area':ccr_area,'ccr_id':ccr_id}   

    
   
# Diagnostics
idiag = False
if idiag:
    for key in ccr_core_dict.keys():
        print(key,', # of times:',len(ccr_core_dict[key].keys()))
        for key2 in ccr_core_dict[key].keys():
            print('Time:',key2,'; # of cores:',len(ccr_core_dict[key][key2]['ccr_id']))
        print(aaaa)
    
print(aaaa)
    
#===============================================
# Construct time series
#===============================================

ccr_core_number = {}
ccr_core_mean_area = {}
ccr_core_tot_area = {}

for key in ccr_core_dict.keys():
    print(key)
    mean_area_arr = []
    tot_area_arr = []
    num_arr = []
    # Loops through time
    for key2 in ccr_core_dict[key].keys():
        tmp_area = ccr_core_dict[key][key2]['ccr_area']
        tmp_id = ccr_core_dict[key][key2]['ccr_id']
        num_arr.append(len(tmp_id))
        mean_area_arr.append(np.mean(tmp_area))
        tot_area_arr.append(np.sum(tmp_area))
    num_arr = np.array(num_arr)
    mean_area_arr = np.array(mean_area_arr)
    tot_area_arr = np.array(tot_area_arr)
    ccr_core_number[key] = num_arr
    ccr_core_mean_area[key] = mean_area_arr
    ccr_core_tot_area[key] = tot_area_arr
    
#=========================================
# perform running mean
#=========================================
ccr_core_number_run_mean_dict = {}
ccr_core_mean_area_run_mean_dict = {}
ccr_core_tot_area_run_mean_dict = {}

irun_mean = 1
if irun_mean == 1:
    N=4
    # Number
    for key,val in ccr_core_number.items():
        tmpval = val.copy()
        tmpval_rm = running_mean(tmpval,N)
        ccr_core_number_run_mean_dict[key] = tmpval_rm
    # Mean Area
    for key,val in ccr_core_mean_area.items():
        tmpval = val.copy()
        tmpval_rm = running_mean(tmpval,N)
        ccr_core_mean_area_run_mean_dict[key] = tmpval_rm    
    # Tot Area
    for key,val in ccr_core_tot_area.items():
        tmpval = val.copy()
        tmpval_rm = running_mean(tmpval,N)
        ccr_core_tot_area_run_mean_dict[key] = tmpval_rm  
        
        
#===============================================
# Make Fig. 8 of manuscript
#===============================================
stoch_short_colors = ['lightsalmon','salmon','red','darkred','firebrick']
stoch_long_colors = ['powderblue','deepskyblue','dodgerblue','blue','navy']
theta_pert_colors = [ 'orchid','fuchsia','mediumorchid','purple','darkorchid'] 
    
#--------------------------------
# Start figure
#--------------------------------

fig = plt.figure(figsize=(12,10))
ax1 = fig.add_subplot(331)
ax2 = fig.add_subplot(332)
ax3 = fig.add_subplot(333)
ax4 = fig.add_subplot(334)
ax5 = fig.add_subplot(335)
ax6 = fig.add_subplot(336)
ax7 = fig.add_subplot(337)
ax8 = fig.add_subplot(338)
ax9 = fig.add_subplot(339)
Fontsize=18

axlist = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]

fac = 1.1
ax1.set_ylabel('Number of CCR\nCores',fontsize=Fontsize*fac)
ax4.set_ylabel('Mean CCR Core\nArea [km$^{2}$]',fontsize=Fontsize*fac)
ax7.set_ylabel('Total CCR Core\nArea [km$^{2}$]',fontsize=Fontsize*fac)
ax7.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize*fac)
ax8.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize*fac)
ax9.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize*fac)

fac2 = 0.9
for ax in axlist:
    #ax.grid(ls='dotted',lw=1,c='grey')
    ax.tick_params('x',labelsize=Fontsize*fac2)
    ax.tick_params('y',labelsize=Fontsize)
    ax.xaxis.set_major_formatter(dfmt)
    #ax.xaxis.set_ticks(time[::48])
    ax.set_xlim(time.min(),time.max())
    ax.xaxis.set_ticks(time[::48])

# remove x labels where we don't need them
dum_axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
for ax in dum_axlist:
    ax.xaxis.set_ticklabels([])
# remove y label where we don't need them
dum_axlist = [ax2,ax3,ax5,ax6,ax8,ax9]
for ax in dum_axlist:
    ax.yaxis.set_ticklabels([])

# Set xticks on bottom row
dum_axlist = [ax7,ax8,ax9]
#vfor ax in dum_axlist:
for ax in axlist:
    ax.xaxis.set_ticks(time[::48])
    
# Add grid lines
for ax in axlist:
    ax.axvline(time[24],lw=1,ls='dotted',c='dimgrey')
    ax.axvline(time[48],lw=1,ls='dotted',c='dimgrey')
    ax.axvline(time[72],lw=1,ls='dotted',c='dimgrey')
    ax.axvline(time[96],lw=1,ls='dotted',c='dimgrey')

# Plot BASEILNE & SPOL on all subplots

# Core Number
axlist1 = [ax1,ax2,ax3]
for ax in axlist1:
    ax.plot(time,ccr_core_number_run_mean_dict['baseline'],c='k',lw=3)
    ax.plot(time,ccr_core_number_run_mean_dict['spol'],c='grey',lw=3)
    ax.set_ylim(0,100)
    ax.axhline(25,lw=1,ls='dotted',c='dimgrey')
    ax.axhline(50,lw=1,ls='dotted',c='dimgrey')
    ax.axhline(75,lw=1,ls='dotted',c='dimgrey')
    
# Add STOCH_SHORT
ax1.plot(time,ccr_core_number_run_mean_dict['stoch1_short'],lw=1.5,c=stoch_short_colors[0])
ax1.plot(time,ccr_core_number_run_mean_dict['stoch2_short'],lw=1.5,c=stoch_short_colors[1])
ax1.plot(time,ccr_core_number_run_mean_dict['stoch3_short'],lw=1.5,c=stoch_short_colors[2])
ax1.plot(time,ccr_core_number_run_mean_dict['stoch4_short'],lw=1.5,c=stoch_short_colors[3])
ax1.plot(time,ccr_core_number_run_mean_dict['stoch5_short'],lw=1.5,c=stoch_short_colors[4])
# Add STOCH_LONG
ax2.plot(time,ccr_core_number_run_mean_dict['stoch1_long'],lw=1.5,c=stoch_long_colors[0])
ax2.plot(time,ccr_core_number_run_mean_dict['stoch2_long'],lw=1.5,c=stoch_long_colors[1])
ax2.plot(time,ccr_core_number_run_mean_dict['stoch3_long'],lw=1.5,c=stoch_long_colors[2])
ax2.plot(time,ccr_core_number_run_mean_dict['stoch4_long'],lw=1.5,c=stoch_long_colors[3])
ax2.plot(time,ccr_core_number_run_mean_dict['stoch5_long'],lw=1.5,c=stoch_long_colors[4])
# THETA_PERT
ax3.plot(time,ccr_core_number_run_mean_dict['theta_pert1'],lw=1.5,c=theta_pert_colors[0])
ax3.plot(time,ccr_core_number_run_mean_dict['theta_pert2'],lw=1.5,c=theta_pert_colors[1])
ax3.plot(time,ccr_core_number_run_mean_dict['theta_pert3'],lw=1.5,c=theta_pert_colors[2])
ax3.plot(time,ccr_core_number_run_mean_dict['theta_pert4'],lw=1.5,c=theta_pert_colors[3])
ax3.plot(time,ccr_core_number_run_mean_dict['theta_pert5'],lw=1.5,c=theta_pert_colors[4])
# Add other sims
ax3.plot(time,ccr_core_number_run_mean_dict['4x'],lw=3,c='k',ls='dashed')
ax3.plot(time,ccr_core_number_run_mean_dict['no_mixing'],lw=3,c='k',ls='dotted')
ax3.plot(time,ccr_core_number_run_mean_dict['baseline_1km_CG'],lw=3,c='darkorange',ls='-.')

# Mean Core Area
axlist2 = [ax4,ax5,ax6]
for ax in axlist2:
    ax.plot(time,ccr_core_mean_area_run_mean_dict['baseline'],c='k',lw=3)
    ax.plot(time,ccr_core_mean_area_run_mean_dict['spol'],c='grey',lw=3)
    ax.set_ylim(0,220)
    ax.axhline(50,lw=1,ls='dotted',c='dimgrey')
    ax.axhline(100,lw=1,ls='dotted',c='dimgrey')
    ax.axhline(150,lw=1,ls='dotted',c='dimgrey')  
    
# Add STOCH_SHORT
ax4.plot(time,ccr_core_mean_area_run_mean_dict['stoch1_short'],lw=1.5,c=stoch_short_colors[0])
ax4.plot(time,ccr_core_mean_area_run_mean_dict['stoch2_short'],lw=1.5,c=stoch_short_colors[1])
ax4.plot(time,ccr_core_mean_area_run_mean_dict['stoch3_short'],lw=1.5,c=stoch_short_colors[2])
ax4.plot(time,ccr_core_mean_area_run_mean_dict['stoch4_short'],lw=1.5,c=stoch_short_colors[3])
ax4.plot(time,ccr_core_mean_area_run_mean_dict['stoch5_short'],lw=1.5,c=stoch_short_colors[4])
# Add STOCH_LONG
ax5.plot(time,ccr_core_mean_area_run_mean_dict['stoch1_long'],lw=1.5,c=stoch_long_colors[0])
ax5.plot(time,ccr_core_mean_area_run_mean_dict['stoch2_long'],lw=1.5,c=stoch_long_colors[1])
ax5.plot(time,ccr_core_mean_area_run_mean_dict['stoch3_long'],lw=1.5,c=stoch_long_colors[2])
ax5.plot(time,ccr_core_mean_area_run_mean_dict['stoch4_long'],lw=1.5,c=stoch_long_colors[3])
ax5.plot(time,ccr_core_mean_area_run_mean_dict['stoch5_long'],lw=1.5,c=stoch_long_colors[4])
# THETA_PERT
ax6.plot(time,ccr_core_mean_area_run_mean_dict['theta_pert1'],lw=1.5,c=theta_pert_colors[0])
ax6.plot(time,ccr_core_mean_area_run_mean_dict['theta_pert2'],lw=1.5,c=theta_pert_colors[1])
ax6.plot(time,ccr_core_mean_area_run_mean_dict['theta_pert3'],lw=1.5,c=theta_pert_colors[2])
ax6.plot(time,ccr_core_mean_area_run_mean_dict['theta_pert4'],lw=1.5,c=theta_pert_colors[3])
ax6.plot(time,ccr_core_mean_area_run_mean_dict['theta_pert5'],lw=1.5,c=theta_pert_colors[4])
# Add other sims
ax6.plot(time,ccr_core_mean_area_run_mean_dict['4x'],lw=3,c='k',ls='dashed')
ax6.plot(time,ccr_core_mean_area_run_mean_dict['no_mixing'],lw=3,c='k',ls='dotted')
ax6.plot(time,ccr_core_mean_area_run_mean_dict['baseline_1km_CG'],lw=3,c='darkorange',ls='-.')    
    
# Total Core Area
axlist3 = [ax7,ax8,ax9]
for ax in axlist3:
    ax.plot(time,ccr_core_tot_area_run_mean_dict['baseline'],c='k',lw=3)
    ax.plot(time,ccr_core_tot_area_run_mean_dict['spol'],c='grey',lw=3)
    ax.set_ylim(0,8000)
    ax.axhline(2000,lw=1,ls='dotted',c='dimgrey')
    ax.axhline(4000,lw=1,ls='dotted',c='dimgrey')
    ax.axhline(6000,lw=1,ls='dotted',c='dimgrey')       
# Add STOCH_SHORT
ax7.plot(time,ccr_core_tot_area_run_mean_dict['stoch1_short'],lw=1.5,c=stoch_short_colors[0])
ax7.plot(time,ccr_core_tot_area_run_mean_dict['stoch2_short'],lw=1.5,c=stoch_short_colors[1])
ax7.plot(time,ccr_core_tot_area_run_mean_dict['stoch3_short'],lw=1.5,c=stoch_short_colors[2])
ax7.plot(time,ccr_core_tot_area_run_mean_dict['stoch4_short'],lw=1.5,c=stoch_short_colors[3])
ax7.plot(time,ccr_core_tot_area_run_mean_dict['stoch5_short'],lw=1.5,c=stoch_short_colors[4])
# Add STOCH_LONG
ax8.plot(time,ccr_core_tot_area_run_mean_dict['stoch1_long'],lw=1.5,c=stoch_long_colors[0])
ax8.plot(time,ccr_core_tot_area_run_mean_dict['stoch2_long'],lw=1.5,c=stoch_long_colors[1])
ax8.plot(time,ccr_core_tot_area_run_mean_dict['stoch3_long'],lw=1.5,c=stoch_long_colors[2])
ax8.plot(time,ccr_core_tot_area_run_mean_dict['stoch4_long'],lw=1.5,c=stoch_long_colors[3])
ax8.plot(time,ccr_core_tot_area_run_mean_dict['stoch5_long'],lw=1.5,c=stoch_long_colors[4])
# THETA_PERT
ax9.plot(time,ccr_core_tot_area_run_mean_dict['theta_pert1'],lw=1.5,c=theta_pert_colors[0])
ax9.plot(time,ccr_core_tot_area_run_mean_dict['theta_pert2'],lw=1.5,c=theta_pert_colors[1])
ax9.plot(time,ccr_core_tot_area_run_mean_dict['theta_pert3'],lw=1.5,c=theta_pert_colors[2])
ax9.plot(time,ccr_core_tot_area_run_mean_dict['theta_pert4'],lw=1.5,c=theta_pert_colors[3])
ax9.plot(time,ccr_core_tot_area_run_mean_dict['theta_pert5'],lw=1.5,c=theta_pert_colors[4])
# Add other sims
ax9.plot(time,ccr_core_tot_area_run_mean_dict['4x'],lw=3,c='k',ls='dashed')
ax9.plot(time,ccr_core_tot_area_run_mean_dict['no_mixing'],lw=3,c='k',ls='dotted')
ax9.plot(time,ccr_core_tot_area_run_mean_dict['baseline_1km_CG'],lw=3,c='darkorange',ls='-.')  



dumx = 0.5
dumy = 1.175
ax1.text(dumx,dumy,'STOCH\_SHORT',fontsize=Fontsize*1.5,fontweight='bold',\
         bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
         rotation=0,va='center',ha='center')
ax2.text(dumx,dumy,'STOCH\_LONG',fontsize=Fontsize*1.5,fontweight='bold',\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
         rotation=0,va='center',ha='center')
ax3.text(dumx,dumy,'$\\theta$-pert',fontsize=Fontsize*1.5,fontweight='bold',\
         bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
         rotation=0,va='center',ha='center')  

custom_lines = [Line2D([0],[0],color = 'black',lw=3),\
                    Line2D([0],[0],color='grey',lw=3),\
                    Line2D([0],[0],color='black',ls='dotted',lw=3),\
                    Line2D([0],[0],color='black',ls='dashed',lw=3),\
                    Line2D([0],[0],color='darkorange',ls='-.',lw=3),\
                   ]

ax8.legend(custom_lines,['BASELINE','SPOL','NO\_MIXING','4X','BASELINE\_1KM\_CG'],fontsize=Fontsize*1.1,\
                    loc='lower center',bbox_to_anchor=(0.5,-0.85),ncol=3,framealpha=True)

labs = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)']
dum = 0
for ax in axlist:
    ax.text(0.24,0.8,labs[dum],transform=ax.transAxes,\
            fontweight='bold',fontsize=Fontsize*2,ha='right')
    dum+=1
    


plt.subplots_adjust(wspace=0.2,hspace=0.2)

outfile = 'fig_08.png'
plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
#plt.show()
plt.close()






