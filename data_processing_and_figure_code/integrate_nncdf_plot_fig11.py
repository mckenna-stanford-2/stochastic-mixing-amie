#==========================================================================
# Title:integrate_nncdf_plot_fig11.py
# Author: McKenna W. Stanford
# Utility: Computes CDFs of nearest neighbor distances (read in by dictinoaries)
# and integrates them. Then computes the organization index (I_org).
# Produces Fig. 11 of manuscript
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
import pandas as pd
from scipy.ndimage import label
import cartopy.crs as ccrs
from matplotlib.cm import get_cmap
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

sim_names = ['baseline',\
             'baseline_1km_CG',\
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

nnd_int_dict = {}

path = '/glade/scratch/mckenna/AMIE/amie_post/updraft_nearest_neighbor/'

# Function to perform CDF
def cdf(data,delta_h):
    data_size=len(data)
    # Set bins edges
    data_set=sorted(set(data))
    if delta_h == 3.:
        bins = np.arange(0,402,2)
    elif delta_h == 1.:
        bins = np.arange(0,202,2)
    # Use the histogram function to bin the data
    counts, bin_edges = np.histogram(data, bins=bins, density=False)
    counts=counts.astype(float)/data_size
    # Find the cdf
    cdf = np.cumsum(counts)
    return cdf, bin_edges

for sim_ii in range(num_sims):
    print('Simulation:',sim_names[sim_ii])
    
    file = path+sim_names[sim_ii]+'_updraft_nearest_neighbor.p'
    pkl_file = open(file,'rb')
    nnd_dict = pickle.load(pkl_file)
    pkl_file.close()
    keys = list(nnd_dict.keys())
    nt = len(keys)
    
    if sim_names[sim_ii] == 'baseline_1km_CG':
        delta_h = 1.
    else:
        delta_h = 3.
    
    nnd_int_dict[sim_names[sim_ii]] = []
    
    for tt in range(nt):
        #print('Time step: {}/{}'.format(tt+1,nt))
        nnd_arr = nnd_dict[str(tt+1)]
        x = nnd_arr*3.
        # Creat CDF
        y,bin_edges = cdf(x,delta_h)
        
        mid_bins = []
        for ii in range(len(bin_edges)-1):
            mid_bins.append(0.5*(bin_edges[ii+1]+bin_edges[ii]))

        r = np.sort(mid_bins)
        if delta_h == 3.:
            L = 400*3.
        else:
            L = 200.*3.
        N = len(x)
        lam = N/(L**2.)
        yran = 1 - np.exp(-1.*lam*np.pi*(r**2.))
        int_nncdf = np.trapz(y,yran)
        
        nnd_int_dict[sim_names[sim_ii]].append(int_nncdf)        
        
        
        
# Forgot to add time to dictionary; read it in from elsewhere
file = '/glade/scratch/mckenna/AMIE/amie_post/inst_rain_rates/baseline_inst_rain_rates.p'
pkl_file = open(file,'rb')
dum_dict = pickle.load(pkl_file)
pkl_file.close()
time = dum_dict['time']
        
for key,val in nnd_int_dict.items():
    print(key,np.shape(val),np.max(val),np.min(val))
        
        
        
#---------------------------------
# Plot Organization Index (I_org)
# and make Fig. 11 of manuscript
#---------------------------------
sims1 = ['baseline','stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
sims2 = ['baseline','stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
sims3 = ['baseline','theta_pert1','theta_pert2','theta_pert3','theta_pert4',\
         'theta_pert5','no_mixing','4x','baseline_1km_CG']

lws1 = [3,1.5,1.5,1.5,1.5,1.5]
lws2 = [3,1.5,1.5,1.5,1.5,1.5]
lws3 = [3,1.5,1.5,1.5,1.5,1.5,3,3,3]

colors1 = ['black','lightsalmon','salmon','red','darkred','firebrick']
colors2 = ['black', 'powderblue','deepskyblue','dodgerblue', 'blue','navy']
colors3 = ['black','orchid','fuchsia','mediumorchid','purple','darkorchid','black','black','darkorange'] 

lss1 = ['solid','solid','solid','solid','solid','solid']
lss2 = ['solid','solid','solid','solid','solid','solid']
lss3 = ['solid','solid','solid','solid','solid','solid','dotted','dashed','-.']


#---------------------------------
# Start Plot
#---------------------------------

fig = plt.figure(figsize=(8.3,14))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
axlist = [ax1,ax2,ax3]
Fontsize=20
fac = 1.1
for ax in axlist:
    ax.grid(which='both',ls='dotted',color='dimgrey',lw=1)
    ax.tick_params(labelsize=Fontsize*fac)
    ax.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize*fac)
    ax.set_ylabel('I$_{org}$',fontsize=Fontsize*fac)
    ax.set_ylim(0.69,0.85)
    ax.set_xlim(time.min(),time.max())

N = 4
dum = 0
for sim in sims1:
    val1 = np.array(nnd_int_dict[sim])
    val1 = running_mean(val1,N)
    ax1.plot(time,val1,color=colors1[dum],ls=lss1[dum],lw=lws1[dum])
    dum+=1
    
dum = 0
for sim in sims2:
    val2 = np.array(nnd_int_dict[sim])
    val2 = running_mean(val2,N)
    ax2.plot(time,val2,color=colors2[dum],ls=lss2[dum],lw=lws2[dum])
    dum+=1
    
dum = 0
for sim in sims3:
    val3 = np.array(nnd_int_dict[sim])
    val3 = running_mean(val3,N)
    ax3.plot(time,val3,color=colors3[dum],ls=lss3[dum],lw=lws3[dum])
    dum+=1
    
labs = ['(a)','(b)','(c)']
dum=0
for ax in axlist:
    ax.xaxis.set_major_formatter(dfmt)
    ax.set_xticks(time[::24])
    ax.text(0.,1.125,labs[dum],fontsize=Fontsize*2,transform=ax.transAxes,fontweight='bold')
    ax.yaxis.set_ticks([0.7,0.75,0.8,0.85])
    dum+=1
    
plt.text(0.5,1.05,'STOCH\_SHORT',fontsize=Fontsize*1.5,fontweight='bold',\
         bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
         ha='center',va='bottom')
plt.text(0.5,1.05,'STOCH\_LONG',fontsize=Fontsize*1.5,fontweight='bold',\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
         ha='center',va='bottom')
plt.text(0.5,1.05,'$\\theta$-pert',fontsize=Fontsize*1.5,fontweight='bold',\
         bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
         ha='center',va='bottom')  

plt.subplots_adjust(hspace=0.55,wspace=0.25)

custom_lines = [Line2D([0],[0],color = 'black',lw=3),\
                Line2D([0],[0],color='black',ls='dashed',lw=3),\
                Line2D([0],[0],color='black',ls='dotted',lw=3),\
                Line2D([0],[0],color='darkorange',ls='-.',lw=3),\
               ]
ax3.legend(custom_lines,['BASELINE','4X','NO\_MIXING','BASELINE\_1KM\_CG'],fontsize=Fontsize,\
                loc='lower center',ncol=2,bbox_to_anchor=(0.5,-0.75))

outfile = 'fig_11.png'
plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
#plt.show()
plt.close()     

       
print(aaaa)



#---------------------------------
# Plots
#---------------------------------
import matplotlib.dates as mdates
dfmt = mdates.DateFormatter('%d-%H')

plt.rc('text',usetex=True)
fig = plt.figure(figsize=(5,11))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
axlist = [ax1,ax2,ax3]
for ax in axlist:
    ax.grid(which='both',ls='dotted',color='grey')
    ax.tick_params(labelsize=17.5)
    ax.set_xlabel('Time [UTC]',fontsize=20)
    ax.set_ylabel('I$_{org}$',fontsize=20)

N = 6
dum = 0
for sim in sims1:
    val1 = np.array(nnd_int_dict[sim])
    val1 = running_mean(val1,N)
    ax1.plot(time,val1,color=colors1[dum],ls=lss1[dum],lw=lws1[dum])
    dum+=1
dum = 0
for sim in sims2:
    val2 = np.array(nnd_int_dict[sim])
    val2 = running_mean(val2,N)
    ax2.plot(time,val2,color=colors2[dum],ls=lss2[dum],lw=lws2[dum])
    dum+=1
dum = 0
for sim in sims3:
    val3 = np.array(nnd_int_dict[sim])
    val3 = running_mean(val3,N)
    ax3.plot(time,val3,color=colors3[dum],ls=lss3[dum],lw=lws3[dum])
    dum+=1
    
tmpval = np.array(nnd_int_1km)
tmpval = running_mean(tmpval,N)
#ax3.plot(time,tmpval,color='darkorange',ls='-.',lw=2)

labs = ['(a)','(b)','(c)']
dum=0
for ax in axlist:
    ax.xaxis.set_major_formatter(dfmt)
    ax.set_xticks([time[0],time[24],time[48],time[72],time[96]])
    ax.text(-0.175,1.125,labs[dum],fontsize=40,transform=ax.transAxes,fontweight='bold')
    dum+=1
    
plt.text(0.5,1.05,'STOCH\_SHORT',fontsize=20,fontweight='bold',\
         bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
         ha='center',va='bottom')
plt.text(0.5,1.05,'STOCH\_LONG',fontsize=20,fontweight='bold',\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
         ha='center',va='bottom')
plt.text(0.5,1.05,'$\\theta$-pert',fontsize=20,fontweight='bold',\
         bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
         ha='center',va='bottom')  

plt.subplots_adjust(hspace=0.55,wspace=0.25)

custom_lines = [Line2D([0],[0],color = 'black',lw=3),\
                Line2D([0],[0],color='black',ls='dashed',lw=3),\
                Line2D([0],[0],color='black',ls='dotted',lw=3),\
                #Line2D([0],[0],color='darkorange',ls='-.',lw=3),\
               ]
#ax3.legend(custom_lines,['BASELINE','4X','NO\_MIXING','BASELINE\_1KM'],fontsize=25,\
ax3.legend(custom_lines,['BASELINE','4X','NO\_MIXING'],fontsize=25,\
                loc='lower center',ncol=1,bbox_to_anchor=(0.5,-1.25))

plt.show()
        

#if False:    
if True:    
    nnd_int_dict = {}
    nnd_int_1km = []

    for key in nnd_dict.keys():
        nnd_int_dict[key] = []

    if True:
        nt = len(time)
        for tt in range(nt):
            print(time[tt])

            for key,val in nnd_dict.items():
                print(key)

                x = val[tt]*3.
                y,bin_edges = cdf(x)
                mid_bins = []
                for ii in range(len(bin_edges)-1):
                    mid_bins.append(0.5*(bin_edges[ii+1]+bin_edges[ii]))

                r = np.sort(mid_bins)
                L = 400*3.
                N = len(x)
                lam = N/(L**2.)
                yran = 1 - np.exp(-1.*lam*np.pi*(r**2.))
                #print(aaaaa)
                int_nncdf = np.trapz(y,yran)
                nnd_int_dict[key].append(int_nncdf)
        #-----------------
        # 1km
        #-----------------
        for tt in range(nt):
            x_1km = nnd_1km[tt]
            y_1km,bin_edges_1km = cdf_1km(x_1km)
            mid_bins_1km = []
            for ii in range(len(bin_edges_1km)-1):
                mid_bins_1km.append(0.5*(bin_edges_1km[ii+1]+bin_edges_1km[ii]))
            r_1km = np.sort(mid_bins_1km)
            L_1km = 600.
            N_1km = len(x_1km)
            lam_1km = N_1km/(L_1km**2.)
            yran_1km = 1 - np.exp(-1.*lam_1km*np.pi*(r_1km**2.))
            int_nncdf_1km = np.trapz(y_1km,yran_1km)
            nnd_int_1km.append(int_nncdf_1km)

    print('done') 