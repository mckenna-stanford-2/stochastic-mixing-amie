#==========================================================================
# Title: plot_pert_correlations_fig14.py
# Author: McKenna W. Stanford
# Utility: Plot correlation between subdomain-mean perturbation  and the 
# of relative difference in stochastic simulations from BASELINE.
# Variables here are volumetric rainfall, convective rain rate, upward
# mass flux, and convective updraft mass flux
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



# Make SPOL range ring
x0 = np.arange(-148.5,151.5,3)
y0 = np.arange(-148.5,151.5,3)
X0,Y0 = np.meshgrid(x0,y0)
radial_dist = np.sqrt(X0**2 + Y0**2)


sim_names = ['baseline',\
         'stoch1_long',\
         'stoch2_long',\
         'stoch3_long',\
         'stoch4_long',\
         'stoch5_long',\
          ]

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


#===============================================
# Volumetric rainfall
#===============================================

sim_names = ['baseline',\
         'stoch1_long',\
         'stoch2_long',\
         'stoch3_long',\
         'stoch4_long',\
         'stoch5_long',\
          ]

vol_rainfall_dict = {}

path = '/glade/scratch/mckenna/AMIE/amie_post/acc_gridscale_precip/'

num_sims = len(sim_names)

for sim_ii in range(num_sims):
    print('Simulation:',sim_names[sim_ii])
    
    # Load file
    pkl_file = open(path+sim_names[sim_ii]+'_acc_gridscale_precip.p','rb')
    tmpfile = pickle.load(pkl_file)
    tmp_rain_rate = tmpfile['rain_rate'][47:,:,:]
    tmp_rain_rate = tmp_rain_rate[:,200-50:200+50,200-50:200+50]
    pkl_file.close()

    # set points outside of 150-km radius to NaN
    dumid = np.where(radial_dist > 150.)
    tmp_rain_rate[:,dumid[0],dumid[1]] = np.nan
    #plt.contourf(tmp_rain_rate[0,:,:])

    domain_sum = np.nansum(tmp_rain_rate,axis=(1,2))
    #domain_cumsum = np.cumsum(domain_sum)
    
    # Volumetric precipitatio (mm km^2)
    #rain_rate_sum = domain_cumsum[-1]*9.
    rain_rate_sum = domain_sum*9.
    vol_rainfall_dict[sim_names[sim_ii]] = rain_rate_sum
    
#===============================================
# Convective Rain Rates
#=============================================== 

sim_names = ['baseline',\
         'stoch1_long',\
         'stoch2_long',\
         'stoch3_long',\
         'stoch4_long',\
         'stoch5_long',\
          ]


conv_rr_dict = {}
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
    con3_mask = con3_mask[:,200-50:200+50,200-50:200+50]
    #dbz = conv_strat_dict['dbz']

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
    dumid = np.where(radial_dist > 150.)
    rain_rate[:,dumid[0],dumid[1]] = np.nan
    con3_mask[:,dumid[0],dumid[1]] = np.nan
    
    #plt.contourf(rain_rate[0,:,:])
    #plt.contourf(con3_mask[0,:,:])
    
    conv_rain_rate_mean = []

    for tt in range(nt):
        tmp_rain_rate = rain_rate[tt,:,:]
        tmp_conv_strat_id = con3_mask[tt,:,:]
        conv_id = np.where(tmp_conv_strat_id == 1)
        if np.size(conv_id) > 0.:
            tmp_conv_rain_rate_mean = np.mean(tmp_rain_rate[conv_id])
        else:
            tmp_conv_rain_rate_mean = np.nan
        conv_rain_rate_mean.append(tmp_conv_rain_rate_mean)
    conv_rain_rate_mean = np.array(conv_rain_rate_mean)    
    conv_rr_dict[sim] = conv_rain_rate_mean


#===============================================
# Mass Fluxes
#===============================================     
 
sim_names = ['baseline',\
     'stoch1_long',\
     'stoch2_long',\
     'stoch3_long',\
     'stoch4_long',\
     'stoch5_long',\
      ]


up_mf_dict = {}
conv_up_mf_dict = {}
path = '/glade/scratch/mckenna/AMIE/amie_post/mass_flux/'  
for sim in sim_names:
    print('Simulation:',sim)
       
    pkl_file = open(path+'{}_int_mass_flux_dict.p'.format(sim),'rb')
    int_mf_dict = pickle.load(pkl_file)
    pkl_file.close()
    up_mf = int_mf_dict['pos_up_mf'][200-50:200+50,200-50:200+50,:]
    conv_up_mf = int_mf_dict['conv_up_mf'][200-50:200+50,200-50:200+50,:]

    time = int_mf_dict['time']
    nt = len(time)
    
    # set points outside of 150-km radius to NaN
    dumid = np.where(radial_dist > 150.)
    up_mf[dumid[0],dumid[1],:] = np.nan
    conv_up_mf[dumid[0],dumid[1],:] = np.nan  
    
    #plt.contourf(up_mf[:,:,0])
    #plt.contourf(conv_up_mf[:,:,0])
    up_mf_mean = []
    conv_up_mf_mean = []

    for tt in range(nt):
        tmp_up_mf = up_mf[:,:,tt]
        tmp_conv_up_mf = conv_up_mf[:,:,tt]
        tmp_up_mf_mean = np.nanmean(tmp_up_mf)

        dumid = np.where(tmp_conv_up_mf > 0.)
        tmp_conv_up_mf_mean = np.nanmean(tmp_conv_up_mf[dumid])
        
        up_mf_mean.append(tmp_up_mf_mean)
        conv_up_mf_mean.append(tmp_conv_up_mf_mean)

        
    up_mf_mean = np.array(up_mf_mean)    
    conv_up_mf_mean = np.array(conv_up_mf_mean) 
    
    up_mf_dict[sim] = up_mf_mean
    conv_up_mf_dict[sim] = conv_up_mf_mean
    

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
    
vol_rainfall_rel_diff = {}
conv_rr_rel_diff = {}
up_mf_rel_diff = {}
conv_up_mf_rel_diff = {}

baseline_vol_rainfall = vol_rainfall_dict['baseline']
baseline_conv_rr = conv_rr_dict['baseline']
baseline_up_mf = up_mf_dict['baseline']
baseline_conv_up_mf = conv_up_mf_dict['baseline']

    
for sim in sim_names:
    print(sim)
    key = sim

    # vol rainfall
    tmpval = vol_rainfall_dict[key]
    tmp_rel_diff = ((tmpval-baseline_vol_rainfall)/baseline_vol_rainfall)*100.
    vol_rainfall_rel_diff[sim] = tmp_rel_diff

    # convective rain rate
    tmpval = conv_rr_dict[key]
    tmp_rel_diff = ((tmpval-baseline_conv_rr)/baseline_conv_rr)*100.
    conv_rr_rel_diff[sim] = tmp_rel_diff        
        
    # Upward mass flux
    tmpval = up_mf_dict[key]
    tmp_rel_diff = ((tmpval-baseline_up_mf)/baseline_up_mf)*100.
    up_mf_rel_diff[sim] = tmp_rel_diff   

    # Convective mass flux
    tmpval = conv_up_mf_dict[key]
    tmp_rel_diff = ((tmpval-baseline_conv_up_mf)/baseline_conv_up_mf)*100.
    conv_up_mf_rel_diff[sim] = tmp_rel_diff   
    
    
    
#=============================================================
#=============================================================
# Plot Fig. 14 - correlations between sub-domain-mean perturbations
# and volumetric rainfall, convective rain rate, upward mass
# flux, and convective updraft mass flux
#=============================================================
#============================================================-


sim_names = ['stoch1_long',\
        'stoch2_long',\
        'stoch3_long',\
        'stoch4_long',\
        'stoch5_long',\
       ]  
    
fig = plt.figure(figsize=(12,18))
ax1 = fig.add_subplot(5,4,1)
ax2 = fig.add_subplot(5,4,5)
ax3 = fig.add_subplot(5,4,9)
ax4 = fig.add_subplot(5,4,13)
ax5 = fig.add_subplot(5,4,17)

ax6 = fig.add_subplot(5,4,2)
ax7 = fig.add_subplot(5,4,6)
ax8 = fig.add_subplot(5,4,10)
ax9 = fig.add_subplot(5,4,14)
ax10 = fig.add_subplot(5,4,18)

ax11 = fig.add_subplot(5,4,3)
ax12 = fig.add_subplot(5,4,7)
ax13 = fig.add_subplot(5,4,11)
ax14 = fig.add_subplot(5,4,15)
ax15 = fig.add_subplot(5,4,19)

ax16 = fig.add_subplot(5,4,4)
ax17 = fig.add_subplot(5,4,8)
ax18 = fig.add_subplot(5,4,12)
ax19 = fig.add_subplot(5,4,16)
ax20 = fig.add_subplot(5,4,20)

Fontsize=18

axlist1 = [ax1,ax2,ax3,ax4,ax5]
axlist2 = [ax6,ax7,ax8,ax9,ax10]
axlist3 = [ax11,ax12,ax13,ax14,ax15]
axlist4 = [ax16,ax17,ax18,ax19,ax20]
labs1 = ['(a)','(b)','(c)','(d)','(e)']
labs2 = ['(f)','(g)','(h)','(i)','(j)']
labs3 = ['(k)','(l)','(m)','(n)','(o)']
labs4 = ['(p)','(q)','(r)','(s)','(t)']
frac_of_time = np.arange(0,len(time),1)/(len(time)-1)

cmap = cm.get_cmap('bone_r',len(frac_of_time))
newcolors = cmap(np.linspace(0,1,len(frac_of_time)))

axlist_master = axlist1 + axlist2 + axlist3 + axlist4
labs_master = labs1 + labs2 + labs3 + labs4

dumi = 0
for ax in axlist_master:
    ax.grid(ls='dotted',lw=1,color='dimgrey')
    ax.tick_params(labelsize=Fontsize)
    ax.set_ylim(-100,100)
    ax.set_xlim(-2.5,2.5)
    ax.text(-0.175,1.15,labs_master[dumi],fontsize=Fontsize*1.7,transform=ax.transAxes)
    dumi+=1
                
dums = 60
#--------------------------------
# Volumetric rainfall
#--------------------------------
for sim_ii in range(len(sim_names)):
    tmp_mean_pert = mean_pert_dict[sim_names[sim_ii]]
    tmp_var_rel_diff = vol_rainfall_rel_diff[sim_names[sim_ii]]
    
    # Scatterplot
    axlist1[sim_ii].scatter(tmp_mean_pert,tmp_var_rel_diff,\
               s=dums,c=newcolors,alpha=0.75)
    
    # Compute coefficient of determination (r^2)
    x = tmp_mean_pert
    y = tmp_var_rel_diff 
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
# Convective Rain Rate
#--------------------------------
for sim_ii in range(len(sim_names)):
    tmp_mean_pert = mean_pert_dict[sim_names[sim_ii]]
    tmp_var_rel_diff = conv_rr_rel_diff[sim_names[sim_ii]]
    
    # Scatterplot
    axlist2[sim_ii].scatter(tmp_mean_pert,tmp_var_rel_diff,\
               s=dums,c=newcolors,alpha=0.75)
    
    # Compute coefficient of determination (r^2)
    x = tmp_mean_pert
    y = tmp_var_rel_diff 
    dumid = np.where(~np.isnan(tmp_var_rel_diff))
    x = x[dumid]
    y = y[dumid]

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
# Upward mass flux
#--------------------------------
for sim_ii in range(len(sim_names)):
    tmp_mean_pert = mean_pert_dict[sim_names[sim_ii]]
    tmp_var_rel_diff = up_mf_rel_diff[sim_names[sim_ii]]
    
    # Scatterplot
    axlist3[sim_ii].scatter(tmp_mean_pert,tmp_var_rel_diff,\
               s=dums,c=newcolors,alpha=0.75)
    
    # Compute coefficient of determination (r^2)
    x = tmp_mean_pert
    y = tmp_var_rel_diff
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
# Convective Mass flux
#--------------------------------
for sim_ii in range(len(sim_names)):
    tmp_mean_pert = mean_pert_dict[sim_names[sim_ii]]
    tmp_var_rel_diff = conv_up_mf_rel_diff[sim_names[sim_ii]]
    
    # Scatterplot
    axlist4[sim_ii].scatter(tmp_mean_pert,tmp_var_rel_diff,\
               s=dums,c=newcolors,alpha=0.75)
    
    # Compute coefficient of determination (r^2)
    x = tmp_mean_pert
    y = tmp_var_rel_diff
    m, b = np.polyfit(x, y, deg=1)
    x_sorted = np.array(sorted(x))
    axlist4[sim_ii].plot(x_sorted, m*x_sorted + b,color='blue',lw=4)

    correlation_matrix = np.corrcoef(x,y)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2

    r2_str = str(np.around(r_squared,2))
    axlist4[sim_ii].text(0.5,1.05,'r$^{2}$ = '+r2_str,transform=axlist4[sim_ii].transAxes,\
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
plt.text(0.5,1.4,'Vol. Rainfall',fontsize=Fontsize*1.8,fontweight='bold',\
         bbox=dict(facecolor='grey', alpha=0.25),transform=ax1.transAxes,\
         rotation=0,ha='center',va='bottom')
plt.text(0.5,1.4,'Conv. RR',fontsize=Fontsize*1.8,fontweight='bold',\
         bbox=dict(facecolor='grey', alpha=0.25),transform=ax6.transAxes,\
         rotation=0,ha='center',va='bottom')
plt.text(0.5,1.4,'Pos. $\\langle$MF$\\rangle$',fontsize=Fontsize*1.8,fontweight='bold',\
         bbox=dict(facecolor='grey', alpha=0.25),transform=ax11.transAxes,\
         rotation=0,ha='center',va='bottom')   
plt.text(0.5,1.4,'Conv. $\\langle$MF$\\rangle$',fontsize=Fontsize*1.8,fontweight='bold',\
         bbox=dict(facecolor='grey', alpha=0.25),transform=ax16.transAxes,\
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

plt.text(1.1,-0.4,'Subdomain-Mean Perturbation ($\\bar{\\gamma}$)',fontsize=Fontsize*2,fontweight='bold',\
        rotation=0,ha='center',va='center',transform=ax10.transAxes)    

plt.subplots_adjust(hspace=0.5,wspace=0.3)

outfile = 'fig_14.png'
plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
#plt.show()
plt.close()     

    
    
    




