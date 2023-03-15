#==========================================================================
# Title: plot_bulk_rain_pert_correlations_conv_strat.py
# Author: McKenna W. Stanford
# Origin Date: 05/26/2020
# Date Modified: 05/26/2020
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
# make function that finds the nearest
# element in array to the given value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

savdir = '/glade/scratch/mckenna/AMIE/amie_post/'

x0 = np.arange(-148.5,151.5,3)
y0 = np.arange(-148.5,151.5,3)
X0,Y0 = np.meshgrid(x0,y0)
radial_dist = np.sqrt(X0**2 + Y0**2)


hours = np.arange(0,24,1)
minutes = np.arange(0,60,15)

wtime = []
for hour in hours:
    for minute in minutes:
        wtime.append(datetime.datetime(2011,12,7,hour,minute))
        
hours2 = np.arange(0,12,1)
minutes2 = np.arange(0,60,15)
for hour in hours2:
    for minute in minutes2:
        wtime.append(datetime.datetime(2011,12,8,hour,minute))


wtime.append(datetime.datetime(2011,12,8,12))


    
x97run = 0
if x97run == 1: 
    sim_names = ['baseline',\
             'stoch1_long',\
             'stoch2_long',\
             'stoch3_long',\
             'stoch4_long',\
             'stoch5_long',\
              ]
    
    mf_int_pos_up_mean_time_series = {}
    mf_int_conv_up_1_mean_time_series = {}


    
    for sim in sim_names:
        mf_int_pos_up_mean_time_series[sim] = []
        mf_int_conv_up_1_mean_time_series[sim] = []

        
     # set points outside of 150-km radius to NaN
    xxtmpid = np.where(radial_dist > 150.)
    pkl_file = open(savdir+'mf_int_spatial_dict','rb')
    tmpfile = pickle.load(pkl_file)
    pkl_file.close()  
    mf_int_pos_up_dict = tmpfile['mf_int_pos_up']
    mf_int_conv_up_1_dict = tmpfile['mf_int_conv_up_1']
        
    for sim in sim_names:
        print(sim)
        key = sim+'_new'
        mf_int_pos_up = mf_int_pos_up_dict[key]
        mf_int_conv_up_1 = mf_int_conv_up_1_dict[key]
        time = wtime[48:]
        nt = len(time)
        
        #print(aaaa)
        for tt in range(nt):            
            tmp_mf_int_pos_up = mf_int_pos_up[200-50:200+50,200-50:200+50,tt]
            tmp_mf_int_conv_up_1 = mf_int_conv_up_1[200-50:200+50,200-50:200+50,tt]
            tmp_mf_int_pos_up[xxtmpid] = np.nan
            tmp_mf_int_conv_up_1[xxtmpid] = np.nan
            
            
            mean_mf_int_pos_up = np.nanmean(tmp_mf_int_pos_up)
            mean_mf_int_conv_up_1 = np.nanmean(tmp_mf_int_conv_up_1[tmp_mf_int_conv_up_1 > 0.])
            
            mf_int_pos_up_mean_time_series[sim].append(mean_mf_int_pos_up)
            mf_int_conv_up_1_mean_time_series[sim].append(mean_mf_int_conv_up_1)
        #print(np.max(mf_int_conv_up_1_mean_time_series[key]))

            

    #wtime = wtime[48:]

#print(aaaaa)

        
#============================================================
# Calculate relative difference in these values from the 
# BASELINE
#============================================================
#if True:
if False:
    sim_names = ['stoch1_long',\
             'stoch2_long',\
             'stoch3_long',\
             'stoch4_long',\
             'stoch5_long']
    
    mf_int_pos_up_rel_diff = {}
    mf_int_conv_up_1_rel_diff = {}

    baseline_mf_int_pos_up = np.array(mf_int_pos_up_mean_time_series['baseline'])
    baseline_mf_int_conv_up_1 = np.array(mf_int_conv_up_1_mean_time_series['baseline'])

    
    for sim in sim_names:
        key = sim
        
        # pos up
        tmpval = np.array(mf_int_pos_up_mean_time_series[key])
        tmp_rel_diff = ((tmpval-baseline_mf_int_pos_up)/baseline_mf_int_pos_up)*100.
        mf_int_pos_up_rel_diff[sim] = tmp_rel_diff

        # conv up 1
        tmpval = np.array(mf_int_conv_up_1_mean_time_series[key])
        tmp_rel_diff = ((tmpval-baseline_mf_int_conv_up_1)/baseline_mf_int_conv_up_1)*100.
        mf_int_conv_up_1_rel_diff[sim] = tmp_rel_diff
        
        #print(np.max(tmp_rel_diff))
        
        
#===============================================
# Calculate domain mean and median multiplicative
# factor
#===============================================
stochrun = 0
if stochrun == 1:
    median_mult_fac_dict = {}
    mean_mult_fac_dict = {}
    mean_pert_dict = {}
    median_pert_dict = {}
    
    sims = ['stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long']
    
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
        mean_pert = np.nanmean(tmp,axis=(1,2))
        median_pert = np.nanmedian(tmp,axis=(1,2))        
        median_mult_fac_dict[sim] = median_mult_fac
        mean_mult_fac_dict[sim] = mean_mult_fac
        median_pert_dict[sim] = median_pert
        mean_pert_dict[sim] = mean_pert








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


    
    
#===========================================
# Plot correlations between domain-mean
# perturbation and relative difference
# in positive mass flux and convective
# updraft mass flux
#===========================================
#if False:
if True:

    plt.rc('text',usetex=True)
   
    sims = ['stoch1_long',\
            'stoch2_long',\
            'stoch3_long',\
            'stoch4_long',\
            'stoch5_long',\
           ]  
    
    
    fig = plt.figure(figsize=(10,20))
    ax1 = fig.add_subplot(5,2,1)
    ax2 = fig.add_subplot(5,2,3)
    ax3 = fig.add_subplot(5,2,5)
    ax4 = fig.add_subplot(5,2,7)
    ax5 = fig.add_subplot(5,2,9)
    ax6 = fig.add_subplot(5,2,2)
    ax7 = fig.add_subplot(5,2,4)
    ax8 = fig.add_subplot(5,2,6)
    ax9 = fig.add_subplot(5,2,8)
    ax10 = fig.add_subplot(5,2,10)
    
    axlist1 = [ax1,ax2,ax3,ax4,ax5]
    axlist2 = [ax6,ax7,ax8,ax9,ax10]
    labs1 = ['(a)','(b)','(c)','(d)','(e)']
    labs2 = ['(f)','(g)','(h)','(i)','(j)']
    frac_of_time = np.arange(0,97,1)/96.
    pink = cm.get_cmap('bone_r',len(frac_of_time))
    newcolors = pink(np.linspace(0,1,len(frac_of_time)))   
    
    dum_mult_fac = np.arange(0.01,100.01,0.1)
    

    #==========================
    # Positive Mass Flux
    #==========================
  
    dum = 0
    for ax in axlist1:
        sim = sims[dum]
        print(sim)
        #ax.set_xlabel('$log_{2}$(Domain Median $F$)',fontsize=15)
        ax.set_xlabel('Domain-mean $\\gamma$',fontsize=25)
        ax.set_ylabel('Rel. Diff. in\nDomain-Mean\nPositive MF [\%]',fontsize=25)
        ax.grid(which='both',ls='solid',lw=2,color='grey')
        ax.tick_params(labelsize=25)
        #ax.set_xscale('log')
        #ax.set_xlim(0.1,10)
        ax.set_ylim(-100,100)
        ax.text(-0.2,1.15,labs1[dum],fontsize=35,transform=ax.transAxes)
        
        #tmp_med_mult_fac = np.log2(median_mult_fac_dict[sim])
        tmp_med_mult_fac = mean_pert_dict[sim]
        #print(aaaa)
        tmp_mf_int_pos_up_rel_diff = mf_int_pos_up_rel_diff[sim]
        ax.scatter(tmp_med_mult_fac,tmp_mf_int_pos_up_rel_diff,\
                   s=75,c=newcolors,alpha=0.75)  
        

        
        x = tmp_med_mult_fac
        y = tmp_mf_int_pos_up_rel_diff 

        #m, b = np.polyfit(x, y, deg=1)
        m, b = np.polyfit(x, y, deg=1)
        #m2, m1, b = np.polyfit(x, y, deg=2)
        x_sorted = np.array(sorted(x))
        ax.plot(x_sorted, m*x_sorted + b,color='blue',lw=5)
        #ax.plot(x_sorted, m2*x_sorted**2. + m1*x_sorted + b,color='blue',lw=5)
        
        #correlation_matrix = np.corrcoef(x,np.log10(y))
        correlation_matrix = np.corrcoef(x,y)
        correlation_xy = correlation_matrix[0,1]
        r_squared = correlation_xy**2

        r2_str = str(np.around(r_squared,2))
        ax.text(0.5,1.05,'r$^{2}$ = '+r2_str,transform=ax.transAxes,\
                fontsize=35,fontweight='bold',va='bottom',ha='center')   
        dum+=1
        
    #==========================
    # Convective Updraft Mass Flux
    #==========================
  
    dum = 0
    for ax in axlist2:
        sim = sims[dum]
        print(sim)
        #ax.set_xlabel('$log_{2}$(Domain Median $F$)',fontsize=15)
        ax.set_xlabel('Domain-mean $\\gamma$',fontsize=25)
        ax.set_ylabel('Rel. Diff. in\nDomain-Mean\nConv. Updraft MF [\%]',fontsize=25)
        ax.grid(which='both',ls='solid',lw=2,color='grey')
        ax.tick_params(labelsize=25)
        #ax.set_xscale('log')
        #ax.set_xlim(0.1,10)
        ax.set_ylim(-100,100)
        ax.text(-0.2,1.15,labs2[dum],fontsize=35,transform=ax.transAxes)
        
        #tmp_med_mult_fac = np.log2(median_mult_fac_dict[sim])
        tmp_med_mult_fac = mean_pert_dict[sim]
        
        tmp_mf_int_conv_up_1_rel_diff = mf_int_conv_up_1_rel_diff[sim]
        ax.scatter(tmp_med_mult_fac,tmp_mf_int_conv_up_1_rel_diff,\
                   s=75,c=newcolors,alpha=0.75)  
        

        
        x = tmp_med_mult_fac
        y = tmp_mf_int_conv_up_1_rel_diff 

        #m, b = np.polyfit(x, y, deg=1)
        m, b = np.polyfit(x, y, deg=1)
        #m2, m1, b = np.polyfit(x, y, deg=2)
        x_sorted = np.array(sorted(x))
        ax.plot(x_sorted, m*x_sorted + b,color='blue',lw=5)
        #ax.plot(x_sorted, m2*x_sorted**2. + m1*x_sorted + b,color='blue',lw=5)
        
        #correlation_matrix = np.corrcoef(x,np.log10(y))
        correlation_matrix = np.corrcoef(x,y)
        correlation_xy = correlation_matrix[0,1]
        r_squared = correlation_xy**2

        r2_str = str(np.around(r_squared,2))
        ax.text(0.5,1.05,'r$^{2}$ = '+r2_str,transform=ax.transAxes,\
                fontsize=35,fontweight='bold',va='bottom',ha='center')   
        dum+=1        
        


    
    bounds = frac_of_time
    bounds2 = [0,0.2,0.4,0.6,0.8,1.]
    norm = matplotlib.colors.BoundaryNorm(bounds,pink.N)
    ax99 = fig.add_axes([0.95,0.125,0.03,0.75])
    cb = matplotlib.colorbar.ColorbarBase(ax99,cmap=pink,norm=norm,\
                                  spacing='uniform',\
                                  ticks=bounds2,\
                                  boundaries=bounds,\
                                  #extend='both',\
                                  orientation='vertical')#,\
                                  #format='%1f')
    ax99.tick_params(labelsize=30)
    ax99.set_ylabel('Fraction of 24-hr Period',fontsize=30)

    plt.text(0.5,1.5,'Positive MF',fontsize=33,fontweight='bold',\
             bbox=dict(facecolor='grey', alpha=0.25),transform=ax1.transAxes,\
             rotation=0,ha='center',va='bottom')
    plt.text(0.5,1.5,'Conv. Updraft MF',fontsize=33,fontweight='bold',\
             bbox=dict(facecolor='grey', alpha=0.25),transform=ax6.transAxes,\
             rotation=0,ha='center',va='bottom')

    plt.subplots_adjust(hspace=0.7,wspace=0.75)
    print(aaaa)       
    
    
    
    
    
