#==========================================================================
# Title: plot_int_mass_flux_fig10.py
# Author: McKenna W. Stanford
# Utility: Plots time series of domain-mean vertically integrated mass
# flux.
# Used to make Fig. 10 of mansuscript.
#==========================================================================


#------------------------------------------
# Imports
#------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import glob
import xarray
import seaborn as sns
import matplotlib
from matplotlib.cm import get_cmap
from netCDF4 import Dataset
import wrf
import pickle
from matplotlib.lines import Line2D
import datetime
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
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
         'baseline_1km_CG',\
          ]

num_sims = len(sim_names)

path = '/glade/scratch/mckenna/AMIE/amie_post/mass_flux/'

# Loop through simulations and compute domain-means

pos_up_mean = {}
conv_up_mean = {}
neg_down_mean = {}
for sim_ii in range(num_sims):
    print('Simulation:',sim_names[sim_ii])
    
    pkl_file = open(path+sim_names[sim_ii]+'_int_mass_flux_dict.p','rb')
    tmpfile = pickle.load(pkl_file)
    pkl_file.close()
    
    pos_up_mean[sim_names[sim_ii]] = np.mean(tmpfile['pos_up_mf'],axis=(0,1))
    neg_down_mean[sim_names[sim_ii]] = np.mean(tmpfile['neg_down_mf'],axis=(0,1))
    conv_up_mf = tmpfile['conv_up_mf']
    nt = len(tmpfile['time'])
    # Loop through times
    conv_up_mean_arr = []
    for tt in range(nt):
        tmp_conv_up_mf = conv_up_mf[:,:,tt]
        dumid = np.where(tmp_conv_up_mf > 0.)
        dum_mean = np.mean(tmp_conv_up_mf[dumid])
        conv_up_mean_arr.append(dum_mean)
    conv_up_mean_arr = np.array(conv_up_mean_arr)
    conv_up_mean[sim_names[sim_ii]] = conv_up_mean_arr
    
    if sim_ii == 0.:
        time = np.array(tmpfile['time'])
    


#-------------------------------------------------
# Plot 
#-------------------------------------------------

sims1 = ['baseline','stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
sims2 = ['baseline','stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
sims3 = ['baseline','theta_pert1','theta_pert2','theta_pert3','theta_pert4',\
         'theta_pert5','no_mixing','4x','baseline_1km_CG']

#lws1 = [3,1.5,1.5,1.5,1.5,1.5]
lws1 = [3,1,1,1,1,1]
#lws2 = [3,1.5,1.5,1.5,1.5,1.5]
lws2 = [3,1,1,1,1,1]
#lws3 = [3,1.5,1.5,1.5,1.5,1.5,3,3]
lws3 = [3,1,1,1,1,1,3,3,3]

colors1 = ['black','lightsalmon','salmon','red','darkred','firebrick']
colors2 = ['black', 'powderblue','deepskyblue','dodgerblue', 'blue','navy']
colors3 = ['black','orchid','fuchsia','mediumorchid','purple','darkorchid','black','black','darkorange'] 

lss1 = ['solid','solid','solid','solid','solid','solid']
lss2 = ['solid','solid','solid','solid','solid','solid']
lss3 = ['solid','solid','solid','solid','solid','solid','dotted','dashed','-.']

#=========================================
# perform running mean
#=========================================
pos_up_mean_run_mean_dict = pos_up_mean.copy()
conv_up_mean_run_mean_dict = conv_up_mean.copy()
neg_down_mean_run_mean_dict = neg_down_mean.copy()

irun_mean = 1
if irun_mean == 1:
    N=4
    # Positive up
    for key,val in pos_up_mean_run_mean_dict.items():
        tmpval = val.copy()
        tmpval_rm = running_mean(tmpval,N)
        pos_up_mean_run_mean_dict[key] = tmpval_rm
    # Negative down  
    for key,val in neg_down_mean_run_mean_dict.items():
        tmpval = val.copy()
        tmpval_rm = running_mean(tmpval,N)
        neg_down_mean_run_mean_dict[key] = tmpval_rm  
    # Conv up
    for key,val in conv_up_mean_run_mean_dict.items():
        tmpval = val.copy()
        tmpval_rm = running_mean(tmpval,N)
        conv_up_mean_run_mean_dict[key] = tmpval_rm  

        
        
        
#-------------------------------------------------
# Start Figure
#-------------------------------------------------
time2 = np.array([pd.to_datetime(time[dd]) for dd in range(len(time))])

fig = plt.figure(figsize=(12,6))
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)
Fontsize=18
labs = ['(a)','(b)','(c)','(d)','(e)','(f)']
axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
dum = 0
for ax in axlist:
    ax.grid(which='both',color='grey',ls='solid',lw=1)
    ax.tick_params(labelsize=Fontsize)
    ax.set_xlim(time[0],time[-1])
    ax.xaxis.set_major_formatter(dfmt)
    ax.text(0.01,0.81,labs[dum],fontsize=Fontsize*1.9,transform=ax.transAxes)
    ax.set_xticks(time2[::48])
    dum+=1
ax1.set_ylabel('Domain-Mean\nUpward $\\langle$MF$\\rangle$\n[x10$^{2}$ kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize)
ax4.set_ylabel('Conditional Mean\nConvective $\\langle$MF$\\rangle$\n[x10$^{2}$ kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize)
ax4.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize)
ax5.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize)
ax6.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize)

axlist1 = [ax1,ax2,ax3]
for ax in axlist1:
    ax.set_ylim(1.8,6.2)
axlist2 = [ax4,ax5,ax6]
for ax in axlist2:
    ax.set_ylim(18,62)

# Upward MF
fac = 1.e-2
dum = 0
for sim in sims1:
    ax1.plot(time,pos_up_mean_run_mean_dict[sim]*fac,color=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    dum+=1
dum = 0
for sim in sims2:
    ax2.plot(time,pos_up_mean_run_mean_dict[sim]*fac,color=colors2[dum],lw=lws2[dum],ls=lss2[dum])
    dum+=1
dum = 0
for sim in sims3:
    ax3.plot(time,pos_up_mean_run_mean_dict[sim]*fac,color=colors3[dum],lw=lws3[dum],ls=lss3[dum])
    dum+=1
    
# Conv. MF
dum = 0
for sim in sims1:
    ax4.plot(time,conv_up_mean_run_mean_dict[sim]*fac,color=colors1[dum],lw=lws1[dum],ls=lss1[dum])
    dum+=1
dum = 0
for sim in sims2:
    ax5.plot(time,conv_up_mean_run_mean_dict[sim]*fac,color=colors2[dum],lw=lws2[dum],ls=lss2[dum])
    dum+=1
dum = 0
for sim in sims3:
    ax6.plot(time,conv_up_mean_run_mean_dict[sim]*fac,color=colors3[dum],lw=lws3[dum],ls=lss3[dum])
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


custom_lines = [Line2D([0],[0],color = 'black',lw=3),\
                Line2D([0],[0],color='black',lw=3,ls='dotted'),\
                Line2D([0],[0],color='black',ls='dashed',lw=3),\
                Line2D([0],[0],color='darkorange',ls='-.',lw=3),\
               ]

ax6.legend(custom_lines,['BASELINE','NO\_MIXING','4X','BASELINE\_1KM\_CG'],fontsize=Fontsize,\
                loc='lower center',bbox_to_anchor=(0.,-0.85),ncol=2)

plt.subplots_adjust(wspace=0.225,hspace=0.25)

outfile = 'fig_10.png'
plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
#plt.show()
plt.close()

print(aaaa)

#----------------------------------------
# w > 2
#----------------------------------------
iiplot = 0
if iiplot == 1:

    N = 4
    fig = plt.figure(figsize=(22,10))
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)
    labs = ['(a)','(b)','(c)','(d)','(e)','(f)']
    axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
    dum = 0
    for ax in axlist:
        ax.grid(which='both',color='grey',ls='solid',lw=1)
        ax.set_xlabel('Time [UTC]',fontsize=25)
        ax.tick_params(labelsize=25)
        #ax.set_ylim(0,17.5)
        ax.set_xlim(time[0],time[-1])
        ax.xaxis.set_major_formatter(dfmt)
        ax.text(0.05,0.825,labs[dum],fontsize=50,transform=ax.transAxes)
        ax.set_xticks([time[0],time[24],time[48],time[72],time[96]])
        dum+=1

    axlist1 = [ax1,ax2,ax3]
    axlist2 = [ax4,ax5,ax6]
    for ax in axlist:
        ax.set_ylabel('Domain-Mean\nVertically-Integrated\nConvective Mass Flux\n[kg m$^{-1}$ s$^{-2}$]',fontsize=25)
    for ax in axlist2:
        ax.set_ylabel('Domain-Mean Conditional\nVertically-Integrated\nConvective Mass Flux\n[kg m$^{-1}$ s$^{-2}$]',fontsize=25)  
        
        
    sims1 = ['baseline','stoch1_short','stoch2_short','stoch3_short','stoch4_short','stoch5_short']
    sims2 = ['baseline','stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
    sims3 = ['baseline','theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5','no_mixing','4x']


    colors1 = ['black',\
                 'lightsalmon',\
                 'salmon',\
                 'red',\
                 'darkred',\
                 'firebrick']


    lws1 = [4,2,2,2,2,2]
    lws2 = [4,2,2,2,2,2]
    lws3 = [4,2,2,2,2,2,4,4]

    colors2 = ['black',\
                 'powderblue',\
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
                 'black',\
                 'black',\
                 ] 

    lss1 = ['solid','solid','solid','solid','solid','solid']
    lss2 = ['solid','solid','solid','solid','solid','solid']
    lss3 = ['solid','solid','solid','solid','solid','solid','dotted','dashed']



    dum = 0
    for sim in sims1:
        val1 = mf_int_conv_up_2_mean_dict[sim]
        val1 = running_mean(val1,N)
        ax1.plot(time,val1,color=colors1[dum],lw=lws1[dum],ls=lss1[dum])
        val1 = mf_int_conv_up_2_mean_conditional_dict[sim]
        val1 = running_mean(val1,N)
        ax4.plot(time,val1,color=colors1[dum],lw=lws1[dum],ls=lss1[dum])
        dum+=1
    
    dum = 0
    for sim in sims2:
        val1 = mf_int_conv_up_2_mean_dict[sim]
        val1 = running_mean(val1,N)
        ax2.plot(time,val1,color=colors2[dum],lw=lws2[dum],ls=lss2[dum])
        val1 = mf_int_conv_up_2_mean_conditional_dict[sim]
        val1 = running_mean(val1,N)
        ax5.plot(time,val1,color=colors2[dum],lw=lws2[dum],ls=lss2[dum])
        dum+=1  
        
    dum = 0
    for sim in sims3:
        val1 = mf_int_conv_up_2_mean_dict[sim]
        val1 = running_mean(val1,N)
        ax3.plot(time,val1,color=colors3[dum],lw=lws3[dum],ls=lss3[dum])
        val1 = mf_int_conv_up_2_mean_conditional_dict[sim]
        val1 = running_mean(val1,N)
        ax6.plot(time,val1,color=colors3[dum],lw=lws3[dum],ls=lss3[dum])
        dum+=1  
        


    plt.text(0.5,1.05,'STOCH\_SHORT',fontsize=33,fontweight='bold',\
             bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
             ha='center',va='bottom')
    plt.text(0.5,1.05,'STOCH\_LONG',fontsize=33,fontweight='bold',\
             bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
             ha='center',va='bottom')
    plt.text(0.5,1.05,'$\\theta$-pert',fontsize=33,fontweight='bold',\
             bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
             ha='center',va='bottom')


    custom_lines = [Line2D([0],[0],color = 'black',lw=4),\
                    Line2D([0],[0],color='black',lw=4,ls='dotted'),\
                    Line2D([0],[0],color='black',ls='dashed',lw=4)]

    ax5.legend(custom_lines,['BASELINE','NO\_MIXING','4X'],fontsize=25,\
                    loc='lower center',bbox_to_anchor=(0.5,-0.6),ncol=3)

    plt.subplots_adjust(wspace=0.7,hspace=0.35)
    plt.show()

    print(aaaaaa)



#----------------------------------------
# Positive MF
#----------------------------------------
iiiplot = 0
if iiiplot == 1:

    N = 4
    fig = plt.figure(figsize=(22,10))
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)
    labs = ['(a)','(b)','(c)','(d)','(e)','(f)']
    axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
    dum = 0
    for ax in axlist:
        ax.grid(which='both',color='grey',ls='solid',lw=1)
        ax.set_xlabel('Time [UTC]',fontsize=25)
        ax.tick_params(labelsize=25)
        #ax.set_ylim(0,17.5)
        ax.set_xlim(time[0],time[-1])
        ax.xaxis.set_major_formatter(dfmt)
        ax.text(0.05,0.825,labs[dum],fontsize=50,transform=ax.transAxes)
        ax.set_xticks([time[0],time[24],time[48],time[72],time[96]])
        dum+=1

    axlist1 = [ax1,ax2,ax3]
    axlist2 = [ax4,ax5,ax6]
    for ax in axlist:
        ax.set_ylabel('Domain-Mean\nVertically-Integrated\nPositive Mass Flux\n[kg m$^{-1}$ s$^{-2}$]',fontsize=25)
    for ax in axlist2:
        ax.set_ylabel('Domain-Mean Conditional\nVertically-Integrated\nPositive Mass Flux\n[kg m$^{-1}$ s$^{-2}$]',fontsize=25)  
        
        
    sims1 = ['baseline_new','stoch1_short_new','stoch2_short_new','stoch3_short_new','stoch4_short_new','stoch5_short_new']
    sims2 = ['baseline_new','stoch1_long_new','stoch2_long_new','stoch3_long_new','stoch4_long_new','stoch5_long_new']
    sims3 = ['baseline_new','theta_pert1_new','theta_pert2_new','theta_pert3_new','theta_pert4_new','theta_pert5_new','no_mixing_new','4x_new']

    colors1 = ['black',\
                 'lightsalmon',\
                 'salmon',\
                 'red',\
                 'darkred',\
                 'firebrick']


    lws1 = [4,2,2,2,2,2]
    lws2 = [4,2,2,2,2,2]
    lws3 = [4,2,2,2,2,2,4,4]

    colors2 = ['black',\
                 'powderblue',\
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
                 'black',\
                 'black',\
                 ] 

    lss1 = ['solid','solid','solid','solid','solid','solid']
    lss2 = ['solid','solid','solid','solid','solid','solid']
    lss3 = ['solid','solid','solid','solid','solid','solid','dotted','dashed']



    dum = 0
    for sim in sims1:
        val1 = mf_int_pos_up_mean_dict[sim]
        val1 = running_mean(val1,N)
        ax1.plot(time,val1,color=colors1[dum],lw=lws1[dum],ls=lss1[dum])
        val1 = mf_int_pos_up_mean_conditional_dict[sim]
        val1 = running_mean(val1,N)
        ax4.plot(time,val1,color=colors1[dum],lw=lws1[dum],ls=lss1[dum])
        dum+=1
    
    dum = 0
    for sim in sims2:
        val1 = mf_int_pos_up_mean_dict[sim]
        val1 = running_mean(val1,N)
        ax2.plot(time,val1,color=colors2[dum],lw=lws2[dum],ls=lss2[dum])
        val1 = mf_int_pos_up_mean_conditional_dict[sim]
        val1 = running_mean(val1,N)
        ax5.plot(time,val1,color=colors2[dum],lw=lws2[dum],ls=lss2[dum])
        dum+=1  
        
    dum = 0
    for sim in sims3:
        val1 = mf_int_pos_up_mean_dict[sim]
        val1 = running_mean(val1,N)
        ax3.plot(time,val1,color=colors3[dum],lw=lws3[dum],ls=lss3[dum])
        val1 = mf_int_pos_up_mean_conditional_dict[sim]
        val1 = running_mean(val1,N)
        ax6.plot(time,val1,color=colors3[dum],lw=lws3[dum],ls=lss3[dum])
        dum+=1  
        


    plt.text(0.5,1.05,'STOCH\_SHORT',fontsize=33,fontweight='bold',\
             bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
             ha='center',va='bottom')
    plt.text(0.5,1.05,'STOCH\_LONG',fontsize=33,fontweight='bold',\
             bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
             ha='center',va='bottom')
    plt.text(0.5,1.05,'$\\theta$-pert',fontsize=33,fontweight='bold',\
             bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
             ha='center',va='bottom')


    custom_lines = [Line2D([0],[0],color = 'black',lw=4),\
                    Line2D([0],[0],color='black',lw=4,ls='dotted'),\
                    Line2D([0],[0],color='black',ls='dashed',lw=4)]

    ax5.legend(custom_lines,['BASELINE','NO\_MIXING','4X'],fontsize=25,\
                    loc='lower center',bbox_to_anchor=(0.5,-0.6),ncol=3)

    plt.subplots_adjust(wspace=0.7,hspace=0.35)
    plt.show()

    print(aaaaaa)
    
    
#----------------------------------------
# Positive MF & Convective conditionally-sampled
# Mean MF
#----------------------------------------
ivplot = 1
if ivplot == 1:

    N = 4
    fig = plt.figure(figsize=(22,10))
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)
    labs = ['(a)','(b)','(c)','(d)','(e)','(f)']
    axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
    dum = 0
    for ax in axlist:
        ax.grid(which='both',color='grey',ls='solid',lw=1)
        ax.set_xlabel('Time [UTC]',fontsize=25)
        ax.tick_params(labelsize=25)
        #ax.set_ylim(0,17.5)
        ax.set_xlim(time[0],time[-1])
        ax.xaxis.set_major_formatter(dfmt)
        ax.text(0.05,0.825,labs[dum],fontsize=50,transform=ax.transAxes)
        ax.set_xticks([time[0],time[24],time[48],time[72],time[96]])
        dum+=1

    axlist1 = [ax1,ax2,ax3]
    axlist2 = [ax4,ax5,ax6]
    for ax in axlist:
        ax.set_ylabel('Domain Mean\nPositive $\\langle$MF$\\rangle$\n[kg m$^{-1}$ s$^{-2}$]',fontsize=25)
    for ax in axlist2:
        ax.set_ylabel('Conditional Mean\nConvective Updraft\n$\\langle$MF$\\rangle$ [kg m$^{-1}$ s$^{-2}$]',fontsize=25)  
        
        
    sims1 = ['baseline_new','stoch1_short_new','stoch2_short_new','stoch3_short_new','stoch4_short_new','stoch5_short_new']
    sims2 = ['baseline_new','stoch1_long_new','stoch2_long_new','stoch3_long_new','stoch4_long_new','stoch5_long_new']
    sims3 = ['baseline_new','theta_pert1_new','theta_pert2_new','theta_pert3_new','theta_pert4_new','theta_pert5_new','no_mixing_new','4x_new']
    
    colors1 = ['black',\
                 'lightsalmon',\
                 'salmon',\
                 'red',\
                 'darkred',\
                 'firebrick']


    lws1 = [4,2,2,2,2,2]
    lws2 = [4,2,2,2,2,2]
    lws3 = [4,2,2,2,2,2,4,4]

    colors2 = ['black',\
                 'powderblue',\
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
                 'black',\
                 'black',\
                 ] 

    lss1 = ['solid','solid','solid','solid','solid','solid']
    lss2 = ['solid','solid','solid','solid','solid','solid']
    lss3 = ['solid','solid','solid','solid','solid','solid','dotted','dashed']



    dum = 0
    for sim in sims1:
        val1 = mf_int_pos_up_mean_dict[sim]
        val1 = running_mean(val1,N)
        ax1.plot(time,val1,color=colors1[dum],lw=lws1[dum],ls=lss1[dum])
        val1 = mf_int_conv_up_1_mean_conditional_dict[sim]
        val1 = running_mean(val1,N)
        ax4.plot(time,val1,color=colors1[dum],lw=lws1[dum],ls=lss1[dum])
        dum+=1
    
    dum = 0
    for sim in sims2:
        val1 = mf_int_pos_up_mean_dict[sim]
        val1 = running_mean(val1,N)
        ax2.plot(time,val1,color=colors2[dum],lw=lws2[dum],ls=lss2[dum])
        val1 = mf_int_conv_up_1_mean_conditional_dict[sim]
        val1 = running_mean(val1,N)
        ax5.plot(time,val1,color=colors2[dum],lw=lws2[dum],ls=lss2[dum])
        dum+=1  
        
    dum = 0
    for sim in sims3:
        val1 = mf_int_pos_up_mean_dict[sim]
        val1 = running_mean(val1,N)
        ax3.plot(time,val1,color=colors3[dum],lw=lws3[dum],ls=lss3[dum])
        val1 = mf_int_conv_up_1_mean_conditional_dict[sim]
        val1 = running_mean(val1,N)
        ax6.plot(time,val1,color=colors3[dum],lw=lws3[dum],ls=lss3[dum])
        dum+=1  
        


    plt.text(0.5,1.05,'STOCH\_SHORT',fontsize=33,fontweight='bold',\
             bbox=dict(facecolor='red', alpha=0.25),transform=ax1.transAxes,\
             ha='center',va='bottom')
    plt.text(0.5,1.05,'STOCH\_LONG',fontsize=33,fontweight='bold',\
             bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
             ha='center',va='bottom')
    plt.text(0.5,1.05,'$\\theta$-pert',fontsize=33,fontweight='bold',\
             bbox=dict(facecolor='purple', alpha=0.25),transform=ax3.transAxes,\
             ha='center',va='bottom')


    custom_lines = [Line2D([0],[0],color = 'black',lw=4),\
                    Line2D([0],[0],color='black',lw=4,ls='dotted'),\
                    Line2D([0],[0],color='black',ls='dashed',lw=4)]

    ax5.legend(custom_lines,['BASELINE','NO\_MIXING','4X'],fontsize=25,\
                    loc='lower center',bbox_to_anchor=(0.5,-0.6),ncol=3)

    plt.subplots_adjust(wspace=0.7,hspace=0.35)
    plt.show()

    print(aaaaaa)    