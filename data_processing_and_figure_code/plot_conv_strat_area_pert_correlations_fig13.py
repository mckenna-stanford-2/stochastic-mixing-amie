#==========================================================================
# Title: plot_conv_strat_area_pert_correlations_fig13.py
# Author: McKenna W. Stanford
# Utility: Plot correlation between domain median/mean multiplicative
# factor of relative change in stochastic simulations from BASELINE.
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
#------------------------------------------
# Parameters
#------------------------------------------
plt.rc('text',usetex=True)
dfmt = mdates.DateFormatter('%d-%H')

path = '/glade/scratch/mckenna/AMIE/amie_post/conv_strat_id/'



sim_names = ['baseline',\
         'stoch1_long',\
         'stoch2_long',\
         'stoch3_long',\
         'stoch4_long',\
         'stoch5_long',\
          ]

    
tot_echo_area_dict = {}
conv_area_dict = {}
strat_area_dict = {}

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
        tmp_con3_mask = con3_mask[tt,:,:]
        tmp_con3_mask = tmp_con3_mask[200-50:200+50,200-50:200+50]
        
        # set points outside of 150-km radius to NaN
        dumid = np.where(radial_dist > 150.) 
        tmp_con3_mask[dumid] = np.nan
        tmp_con3_mask = np.ndarray.flatten(tmp_con3_mask)
        tmp_con3_mask = tmp_con3_mask[~np.isnan(tmp_con3_mask)]
        
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
    
time = tmpfile['time']
    
# Plot time series
iplot = False
if iplot:
    #time2 = np.array([pd.to_datetime(time[dd]) for dd in range(len(time))])
    #time = time2
    sims = ['baseline','stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
    lws = [3,1,1,1,1,1]
    colors = ['black', 'powderblue','deepskyblue','dodgerblue', 'blue','navy']
    lss = ['solid','solid','solid','solid','solid','solid']
    
    
    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    Fontsize=14
    axlist = [ax1,ax2,ax3]
    
    for ax in axlist:
        ax.set_xlabel('Time',fontsize=Fontsize)
        ax.tick_params(labelsize=Fontsize)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_xticks(time[::48])
        ax.grid()
    ax1.set_ylabel('Total Echo Area [km$^{2}$]',fontsize=Fontsize)
    ax2.set_ylabel('Convective Area [km$^{2}$]',fontsize=Fontsize)
    ax3.set_ylabel('Stratiform Area [km$^{2}$]',fontsize=Fontsize)

    dum = 0
    for sim in sims:
        ax1.plot(time,tot_echo_area_dict[sim],c=colors[dum],lw=lws[dum],ls=lss[dum])
        ax2.plot(time,conv_area_dict[sim],c=colors[dum],lw=lws[dum],ls=lss[dum])
        ax3.plot(time,strat_area_dict[sim],c=colors[dum],lw=lws[dum],ls=lss[dum])
        dum+=1

        
    plt.subplots_adjust()
    plt.show()
    plt.close()    
    
    
        
#============================================================
# Calculate relative difference in these values from the 
# BASELINE
#============================================================

sim_names = ['stoch1_long',\
         'stoch2_long',\
         'stoch3_long',\
         'stoch4_long',\
         'stoch5_long']
    
conv_area_rel_diff = {}
strat_area_rel_diff = {}
tot_echo_area_rel_diff = {}

baseline_conv_area = conv_area_dict['baseline']
baseline_strat_area = strat_area_dict['baseline']
baseline_tot_echo_area = tot_echo_area_dict['baseline']

    
for sim in sim_names:
    print(sim)
    key = sim

    # conv area
    tmpval = conv_area_dict[key]
    tmp_rel_diff = ((tmpval-baseline_conv_area)/baseline_conv_area)*100.
    conv_area_rel_diff[sim] = tmp_rel_diff

     # strat area
    tmpval = np.array(strat_area_dict[key])
    tmp_rel_diff = ((tmpval-baseline_strat_area)/baseline_strat_area)*100.
    strat_area_rel_diff[sim] = tmp_rel_diff

     # tot echo area
    tmpval = np.array(tot_echo_area_dict[key])
    tmp_rel_diff = ((tmpval-baseline_tot_echo_area)/baseline_tot_echo_area)*100.
    tot_echo_area_rel_diff[sim] = tmp_rel_diff      

        
        
#===============================================
# Calculate domain mean and median multiplicative
# factor
#===============================================

icalc = False
if icalc:
    mean_pert_dict = {}

    sim_names = [ 'stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
    num_sims = len(sim_names)
    for sim_ii in range(num_sims):
        print('Simulation:',sim_names[sim_ii])
        tmppath = '/glade/scratch/mckenna/AMIE/'
        files = sorted(glob.glob(tmppath+sim_names[sim_ii]+'/wrfout*'))
        # Limit to after 12 UTC
        files = files[48:]
        num_files = len(files)
        mean_perts = []
        for tt in range(num_files):
            ncfile = xarray.open_dataset(files[tt])
            rand_pert = ncfile['RAND_PERT'].values.squeeze().T
            ncfile.close()
            rand_pert = rand_pert[:,:,0]
            rand_pert = rand_pert[200-50:200+50,200-50:200+50]
            tmpid = np.where(radial_dist > 150.)
            rand_pert[tmpid] = np.nan
            # plt.contourf(rand_pert)
            rand_pert = np.ndarray.flatten(rand_pert)
            rand_pert = rand_pert[~np.isnan(rand_pert)]
            mean_pert = np.nanmean(rand_pert)
            mean_perts.append(mean_pert)
        mean_pert_dict[sim_names[sim_ii]] = np.array(mean_perts)


    # Save to dictionary
    savdir = '/glade/scratch/mckenna/AMIE/amie_post/stats/'
    f = open(savdir+'stoch_long_mean_rand_perts.p','wb')
    pickle.dump(mean_pert_dict,f)
    f.close()
else:
    # Read in dictionary
    tmppath = '/glade/scratch/mckenna/AMIE/amie_post/stats/'
    file = tmppath+'stoch_long_mean_rand_perts.p'
    pkl_file = open(file,'rb')
    mean_pert_dict = pickle.load(pkl_file)




#=============================================================
#=============================================================
# Plot Fig. 13 - correlations between total echo area,
# convective area, and stratiform area and the sub-domain
# mean perturbation
#=============================================================
#============================================================-


sim_names = ['stoch1_long',\
        'stoch2_long',\
        'stoch3_long',\
        'stoch4_long',\
        'stoch5_long',\
       ]  
    
fig = plt.figure(figsize=(12,18))
ax1 = fig.add_subplot(5,3,1)
ax2 = fig.add_subplot(5,3,4)
ax3 = fig.add_subplot(5,3,7)
ax4 = fig.add_subplot(5,3,10)
ax5 = fig.add_subplot(5,3,13)
ax6 = fig.add_subplot(5,3,2)
ax7 = fig.add_subplot(5,3,5)
ax8 = fig.add_subplot(5,3,8)
ax9 = fig.add_subplot(5,3,11)
ax10 = fig.add_subplot(5,3,14)
ax11 = fig.add_subplot(5,3,3)
ax12 = fig.add_subplot(5,3,6)
ax13 = fig.add_subplot(5,3,9)
ax14 = fig.add_subplot(5,3,12)
ax15 = fig.add_subplot(5,3,15)

Fontsize=20

axlist1 = [ax1,ax2,ax3,ax4,ax5]
axlist2 = [ax6,ax7,ax8,ax9,ax10]
axlist3 = [ax11,ax12,ax13,ax14,ax15]
labs1 = ['(a)','(b)','(c)','(d)','(e)']
labs2 = ['(f)','(g)','(h)','(i)','(j)']
labs3 = ['(k)','(l)','(m)','(n)','(o)']
frac_of_time = np.arange(0,len(time),1)/(len(time)-1)

cmap = cm.get_cmap('bone_r',len(frac_of_time))
newcolors = cmap(np.linspace(0,1,len(frac_of_time)))

axlist_master = axlist1 + axlist2 + axlist3
labs_master = labs1 + labs2 + labs3

dumi = 0
for ax in axlist_master:
    ax.grid(ls='dotted',lw=1,color='dimgrey')
    ax.tick_params(labelsize=Fontsize)
    ax.set_ylim(-100,100)
    ax.set_xlim(-2.5,2.5)
    ax.text(-0.175,1.15,labs_master[dumi],fontsize=Fontsize*1.7,transform=ax.transAxes)
    dumi+=1
                
dums=60
    
#--------------------------------
# Total Echo Area
#--------------------------------
for sim_ii in range(len(sim_names)):
    tmp_mean_pert = mean_pert_dict[sim_names[sim_ii]]
    tmp_tot_echo_area_rel_diff = tot_echo_area_rel_diff[sim_names[sim_ii]]
    
    # Scatterplot
    axlist1[sim_ii].scatter(tmp_mean_pert,tmp_tot_echo_area_rel_diff,\
               s=dums,c=newcolors,alpha=0.75)
    
    
    # Compute coefficient of determination (r^2)
    x = tmp_mean_pert
    y = tmp_tot_echo_area_rel_diff 
    m, b = np.polyfit(x, y, deg=1)
    x_sorted = np.array(sorted(x))
    axlist1[sim_ii].plot(x_sorted, m*x_sorted + b,color='blue',lw=4)

    correlation_matrix = np.corrcoef(x,y)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2

    r2_str = str(np.around(r_squared,2))
    axlist1[sim_ii].text(0.5,1.05,'r$^{2}$ = '+r2_str,transform=axlist1[sim_ii].transAxes,\
            fontsize=Fontsize*1.6,fontweight='bold',va='bottom',ha='center',color='blue')   
        
#--------------------------------
# Convective Area
#--------------------------------
for sim_ii in range(len(sim_names)):
    tmp_mean_pert = mean_pert_dict[sim_names[sim_ii]]
    tmp_conv_area_rel_diff = conv_area_rel_diff[sim_names[sim_ii]]
    
    # Scatterplot
    axlist2[sim_ii].scatter(tmp_mean_pert,tmp_conv_area_rel_diff,\
               s=dums,c=newcolors,alpha=0.75)          

    # Compute coefficient of determination (r^2)
    x = tmp_mean_pert
    y = tmp_conv_area_rel_diff 
    m, b = np.polyfit(x, y, deg=1)
    x_sorted = np.array(sorted(x))
    axlist2[sim_ii].plot(x_sorted, m*x_sorted + b,color='blue',lw=4)

    correlation_matrix = np.corrcoef(x,y)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2

    r2_str = str(np.around(r_squared,2))
    axlist2[sim_ii].text(0.5,1.05,'r$^{2}$ = '+r2_str,transform=axlist2[sim_ii].transAxes,\
            fontsize=Fontsize*1.6,fontweight='bold',va='bottom',ha='center',color='blue')  
    
#--------------------------------
# Stratiform Area
#--------------------------------
for sim_ii in range(len(sim_names)):
    tmp_mean_pert = mean_pert_dict[sim_names[sim_ii]]
    tmp_strat_area_rel_diff = strat_area_rel_diff[sim_names[sim_ii]]
    
    # Scatterplot
    axlist3[sim_ii].scatter(tmp_mean_pert,tmp_strat_area_rel_diff,\
               s=dums,c=newcolors,alpha=0.75)           
        

    # Compute coefficient of determination (r^2)
    x = tmp_mean_pert
    y = tmp_strat_area_rel_diff 
    m, b = np.polyfit(x, y, deg=1)
    x_sorted = np.array(sorted(x))
    axlist3[sim_ii].plot(x_sorted, m*x_sorted + b,color='blue',lw=4)

    correlation_matrix = np.corrcoef(x,y)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2

    r2_str = str(np.around(r_squared,2))
    axlist3[sim_ii].text(0.5,1.05,'r$^{2}$ = '+r2_str,transform=axlist3[sim_ii].transAxes,\
            fontsize=Fontsize*1.6,fontweight='bold',va='bottom',ha='center',color='blue')  
      
#--------------------------------
# Plot Color Bar
#--------------------------------
bounds = frac_of_time
dum_ticks = [0,0.2,0.4,0.6,0.8,1.]
norm = matplotlib.colors.BoundaryNorm(bounds,cmap.N)
cbar_ax = fig.add_axes([0.92,0.125,0.03,0.76])
cb = matplotlib.colorbar.ColorbarBase(cbar_ax,cmap=cmap,norm=norm,\
                              spacing='uniform',\
                              ticks=dum_ticks,\
                              boundaries=bounds,\
                              orientation='vertical')
cbar_ax.tick_params(labelsize=Fontsize*1.5)
cbar_ax.set_ylabel('Fraction of 24-hr Period',fontsize=Fontsize*1.5)

#--------------------------------
# Titles
#--------------------------------
plt.text(0.5,1.4,'Total Echo Area',fontsize=Fontsize*1.8,fontweight='bold',\
         bbox=dict(facecolor='grey', alpha=0.25),transform=ax1.transAxes,\
         rotation=0,ha='center',va='bottom')
plt.text(0.5,1.4,'Conv. Area',fontsize=Fontsize*1.8,fontweight='bold',\
         bbox=dict(facecolor='grey', alpha=0.25),transform=ax6.transAxes,\
         rotation=0,ha='center',va='bottom')
plt.text(0.5,1.4,'Strat. Area',fontsize=Fontsize*1.8,fontweight='bold',\
         bbox=dict(facecolor='grey', alpha=0.25),transform=ax11.transAxes,\
         rotation=0,ha='center',va='bottom')   

plt.text(-0.45,0.5,'STOCH1\_LONG',fontsize=Fontsize*1.35,fontweight='bold',\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax1.transAxes,\
         rotation=90,ha='center',va='center')    

plt.text(-0.45,0.5,'STOCH2\_LONG',fontsize=Fontsize*1.35,fontweight='bold',\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax2.transAxes,\
         rotation=90,ha='center',va='center')        

plt.text(-0.45,0.5,'STOCH3\_LONG',fontsize=Fontsize*1.35,fontweight='bold',\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax3.transAxes,\
         rotation=90,ha='center',va='center')   

plt.text(-0.45,0.5,'STOCH4\_LONG',fontsize=Fontsize*1.35,fontweight='bold',\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax4.transAxes,\
         rotation=90,ha='center',va='center')   

plt.text(-0.45,0.5,'STOCH5\_LONG',fontsize=Fontsize*1.35,fontweight='bold',\
         bbox=dict(facecolor='blue', alpha=0.25),transform=ax5.transAxes,\
         rotation=90,ha='center',va='center')    

plt.text(-0.75,0.5,'Relative Difference from BASELINE [\%]',fontsize=Fontsize*2,fontweight='bold',\
        rotation=90,ha='center',va='center',transform=ax3.transAxes)

plt.text(0.5,-0.4,'Subdomain-Mean Perturbation ($\\bar{\\gamma}$)',fontsize=Fontsize*2,fontweight='bold',\
        rotation=0,ha='center',va='center',transform=ax10.transAxes)    

plt.subplots_adjust(hspace=0.4,wspace=0.3)

outfile = 'fig_13.png'
plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
#plt.show()
plt.close()     

   
    
    
    




