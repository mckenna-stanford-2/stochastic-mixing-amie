#==========================================================================
# Title: plot_updraft_tracer_fig12.py
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
from scipy.ndimage import label, find_objects
import cartopy.crs as ccrs
from matplotlib.cm import get_cmap
import pandas as pd
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
#------------------------------------------
# Parameters
#------------------------------------------
plt.rc('text',usetex=True)

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
path = '/glade/scratch/mckenna/AMIE/amie_post/updraft_tracer/'
mean_tracer_dict = {}
mean_w_dict = {}
mean_twc_dict = {}

for sim_ii in range(num_sims):
    print('Simulation:',sim_names[sim_ii])
    
    pkl_file = open(path+sim_names[sim_ii]+'_updraft_tracer.p','rb')
    tmpfile = pickle.load(pkl_file)
    pkl_file.close()  
    

    # Save to diciontary
    mean_tracer_dict[sim_names[sim_ii]] = tmpfile['mean_updraft_tracer']
    mean_w_dict[sim_names[sim_ii]] = tmpfile['mean_updraft_w']
    mean_twc_dict[sim_names[sim_ii]] = tmpfile['mean_updraft_twc']
    z = tmpfile['z']



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

#===========================================
# Plot profiles of avg tracer
#===========================================


fig = plt.figure(figsize=(12,5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
Fontsize=20
axlist = [ax1,ax2,ax3]
for ax in axlist:
    ax.tick_params(labelsize=Fontsize)
    ax.grid(ls='dotted',c='dimgrey',lw=1)
    ax.set_ylabel('Height [km]',fontsize=Fontsize)
    ax.set_ylim(0,17.5)
    ax.set_xlabel('Avg. $\\phi_{updraft}$ [kg$^{-1}$]',fontsize=Fontsize)
    ax.set_xlim(0,0.6)

dum=0
for sim in sims1:
    ax1.plot(mean_tracer_dict[sim],z/1.e3,c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    dum+=1

dum=0
for sim in sims2:
    ax2.plot(mean_tracer_dict[sim],z/1.e3,c=colors2[dum],lw=lws2[dum],ls=lss2[dum])
    dum+=1

dum=0
for sim in sims3:
    ax3.plot(mean_tracer_dict[sim],z/1.e3,c=colors3[dum],lw=lws3[dum],ls=lss3[dum])
    dum+=1

custom_lines = [Line2D([0],[0],color = 'black',lw=3),\
                Line2D([0],[0],color='black',ls='dashed',lw=3),\
                Line2D([0],[0],color='black',ls='dotted',lw=3),\
               ]
ax3.legend(custom_lines,['BASELINE','4X','NO\_MIXING',],fontsize=Fontsize,\
                loc='lower center',ncol=3,bbox_to_anchor=(-0.25,-0.4))
    
plt.text(0.5,1.075,'STOCH\_SHORT',fontsize=Fontsize*1.5,\
         bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
         rotation=0,ha='center')
plt.text(0.5,1.075,'STOCH\_LONG',fontsize=Fontsize*1.5,\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
         rotation=0,ha='center')
plt.text(0.5,1.075,'$\\theta$-pert',fontsize=Fontsize*1.5,\
         bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
         rotation=0,ha='center')

dumx = 0.73
dumy = 0.87
ax1.text(dumx,dumy,'(a)',fontsize=Fontsize*2.,transform=ax1.transAxes)
ax2.text(dumx,dumy,'(b)',fontsize=Fontsize*2.,transform=ax2.transAxes)
ax3.text(dumx,dumy,'(c)',fontsize=Fontsize*2.,transform=ax3.transAxes)

plt.subplots_adjust(wspace=0.35)

outfile = 'fig_12.png'
plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
#plt.show()
plt.close() 


print(aaaa)


#===========================================
# Plot profiles of avg w
#===========================================


fig = plt.figure(figsize=(12,5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
Fontsize=20
axlist = [ax1,ax2,ax3]
for ax in axlist:
    ax.tick_params(labelsize=Fontsize)
    ax.grid(ls='dotted',c='dimgrey',lw=1)
    ax.set_ylabel('Height [km]',fontsize=Fontsize)
    ax.set_ylim(0,17.5)
    ax.set_xlabel('Avg. $w_{updraft}$ [m s$^{-1}$]',fontsize=Fontsize)

dum=0
for sim in sims1:
    ax1.plot(mean_w_dict[sim],z/1.e3,c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    dum+=1

dum=0
for sim in sims2:
    ax2.plot(mean_w_dict[sim],z/1.e3,c=colors2[dum],lw=lws2[dum],ls=lss2[dum])
    dum+=1

dum=0
for sim in sims3:
    ax3.plot(mean_w_dict[sim],z/1.e3,c=colors3[dum],lw=lws3[dum],ls=lss3[dum])
    dum+=1

custom_lines = [Line2D([0],[0],color = 'black',lw=3),\
                Line2D([0],[0],color='black',ls='dashed',lw=3),\
                Line2D([0],[0],color='black',ls='dotted',lw=3),\
               ]
ax3.legend(custom_lines,['BASELINE','4X','NO\_MIXING',],fontsize=Fontsize,\
                loc='lower center',ncol=3,bbox_to_anchor=(-0.25,-0.4))
    
plt.text(0.5,1.075,'STOCH\_SHORT',fontsize=Fontsize*1.5,\
         bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
         rotation=0,ha='center')
plt.text(0.5,1.075,'STOCH\_LONG',fontsize=Fontsize*1.5,\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
         rotation=0,ha='center')
plt.text(0.5,1.075,'$\\theta$-pert',fontsize=Fontsize*1.5,\
         bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
         rotation=0,ha='center')

dumx = 0.73
dumy = 0.87
ax1.text(dumx,dumy,'(a)',fontsize=Fontsize*2.,transform=ax1.transAxes)
ax2.text(dumx,dumy,'(b)',fontsize=Fontsize*2.,transform=ax2.transAxes)
ax3.text(dumx,dumy,'(c)',fontsize=Fontsize*2.,transform=ax3.transAxes)

plt.subplots_adjust(wspace=0.35)

outfile = 'w_profs.png'
#plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
plt.show()
plt.close()     


#===========================================
# Plot profiles of avg TWC
#===========================================


fig = plt.figure(figsize=(12,5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
Fontsize=20
axlist = [ax1,ax2,ax3]
for ax in axlist:
    ax.tick_params(labelsize=Fontsize)
    ax.grid(ls='dotted',c='dimgrey',lw=1)
    ax.set_ylabel('Height [km]',fontsize=Fontsize)
    ax.set_ylim(0,17.5)
    ax.set_xlabel('Avg. TWC$_{updraft}$ [g m$^{-3}$]',fontsize=Fontsize)

dum=0
for sim in sims1:
    ax1.plot(mean_twc_dict[sim],z/1.e3,c=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    dum+=1

dum=0
for sim in sims2:
    ax2.plot(mean_twc_dict[sim],z/1.e3,c=colors2[dum],lw=lws2[dum],ls=lss2[dum])
    dum+=1

dum=0
for sim in sims3:
    ax3.plot(mean_twc_dict[sim],z/1.e3,c=colors3[dum],lw=lws3[dum],ls=lss3[dum])
    dum+=1

custom_lines = [Line2D([0],[0],color = 'black',lw=3),\
                Line2D([0],[0],color='black',ls='dashed',lw=3),\
                Line2D([0],[0],color='black',ls='dotted',lw=3),\
               ]
ax3.legend(custom_lines,['BASELINE','4X','NO\_MIXING',],fontsize=Fontsize,\
                loc='lower center',ncol=3,bbox_to_anchor=(-0.25,-0.4))
    
plt.text(0.5,1.075,'STOCH\_SHORT',fontsize=Fontsize*1.5,\
         bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
         rotation=0,ha='center')
plt.text(0.5,1.075,'STOCH\_LONG',fontsize=Fontsize*1.5,\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
         rotation=0,ha='center')
plt.text(0.5,1.075,'$\\theta$-pert',fontsize=Fontsize*1.5,\
         bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
         rotation=0,ha='center')

dumx = 0.73
dumy = 0.87
ax1.text(dumx,dumy,'(a)',fontsize=Fontsize*2.,transform=ax1.transAxes)
ax2.text(dumx,dumy,'(b)',fontsize=Fontsize*2.,transform=ax2.transAxes)
ax3.text(dumx,dumy,'(c)',fontsize=Fontsize*2.,transform=ax3.transAxes)

plt.subplots_adjust(wspace=0.35)

outfile = 'w_profs.png'
#plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
plt.show()
plt.close()     