#==========================================================================
# Title: compute_updraft_tracer.py
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
z_interps = np.arange(1,17.25,0.25)*1.e3

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

#sim_names = sim_names[0:4]
#sim_names = sim_names[4:8]
#sim_names = sim_names[8:12]
#sim_names = sim_names[12:15]
sim_names = sim_names[15:]


num_sims = len(sim_names)
path = '/glade/scratch/mckenna/AMIE/'
#mean_tracer_dict = {}
#avg_z_dict = {}
#mean_w_dict = {}
#mean_twc_dict = {}

for sim_ii in range(num_sims):
    print('Simulation:',sim_names[sim_ii])
    # Start at 15Z, 3 hrs after tracers are released
    files = sorted(glob.glob(path+sim_names[sim_ii]+'/wrfout*'))[60:]
    nt = len(files)
    
    tracer_all = []
    up_mask_all = []
    time_all = []
    z_all = []
    w_all = []
    twc_all = []
    #========================================================
    # Loop through times
    #========================================================
    #for tt in range(-10,-1):
    for tt in range(nt):
        print('Time Step: {}/{}'.format(tt+1,nt))
        tmpfile = files[tt]

        ncfile = Dataset(tmpfile)
        w = wrf.getvar(ncfile,'wa',meta=False)
        lon = wrf.getvar(ncfile,'lon',meta=False)
        lat = wrf.getvar(ncfile,'lat',meta=False)
        qi = wrf.getvar(ncfile,'QICE',meta=False)
        qr = wrf.getvar(ncfile,'QRAIN',meta=False)
        qc = wrf.getvar(ncfile,'QCLOUD',meta=False)
        tv = wrf.getvar(ncfile,'tv',units='K',meta=False)
        pres = wrf.getvar(ncfile,'pres',units='hPa',meta=False)
        z = wrf.getvar(ncfile,'height',meta=False)
        time = pd.to_datetime(wrf.getvar(ncfile,'times',meta=False))
        tracer = wrf.getvar(ncfile,'TR1',meta=False)
        ncfile.close()
        qt = qi+qr+qc
        rd = 287.04
        rho_air = 100.*pres/(rd*tv)
        twc = rho_air*qt*1.e3
        
        
        i_interp = True
        if i_interp:
            nz_interp = len(z_interps)
            nx = len(twc[0,0,:])
            ny = len(twc[0,:,0])
            
            tracer_interp = np.zeros((nz_interp,ny,nx))
            w_interp = np.zeros((nz_interp,ny,nx))
            twc_interp = np.zeros((nz_interp,ny,nx))
                
            for kk in range(nz_interp):
                dum_tracer = wrf.interplevel(tracer,z,z_interps[kk])
                dum_w = wrf.interplevel(w,z,z_interps[kk])
                dum_twc = wrf.interplevel(twc,z,z_interps[kk])
                tracer_interp[kk,:,:] = dum_tracer
                w_interp[kk,:,:] = dum_w
                twc_interp[kk,:,:] = dum_twc
                
            w = w_interp
            twc = twc_interp
            tracer = tracer_interp


        
        nx = len(w[0,0,:])
        ny = len(w[0,:,0])
        nz = len(w[:,0,0])
        
        # Identify cloudy updrafts
        up_mask = np.zeros(np.shape(w))
        up_id = np.where( (w > 1.) & (twc > 0.1) )
        up_mask[up_id] = 1.
        
        up_mask_all.append(up_mask)
        tracer_all.append(tracer)
        time_all.append(time)
        w_all.append(w)
        twc_all.append(twc)
        
        
    time_all = np.array(time_all)
    tracer_all = np.array(tracer_all)
    up_mask_all = np.array(up_mask_all)
    w_all = np.array(w_all)
    twc_all = np.array(twc_all)
    
    tracer_all = tracer_all.T
    w_all = w_all.T
    twc_all = twc_all.T
    up_mask_all = up_mask_all.T

    nz = len(z_interps)
    nt = len(time_all)
    nx = len(tracer_all[0,:,0,0])
    ny = len(tracer_all[0,0,:,0])
    
    

    # Average across domain and then across time
    if False:
        mean_updraft_tracer_time = np.zeros((nt,nz))
        mean_updraft_w_time = np.zeros((nt,nz))
        mean_updraft_twc_time = np.zeros((nt,nz))
        for tt in range(nt):
            for kk in range(nz):
                #if z_interps[kk] < 1000.:
                #    continue
                tmp_tracer = np.ndarray.flatten(tracer_all[tt,:,:,kk])
                tmp_w = np.ndarray.flatten(w_all[tt,:,:,kk])
                tmp_twc = np.ndarray.flatten(twc_all[tt,:,:,kk])
                tmp_mask = np.ndarray.flatten(up_mask_all[tt,:,:,kk])
                dumid = np.squeeze(np.where(tmp_mask == 1.))
                mean_updraft_tracer_time[tt,kk] = np.mean(tmp_tracer[dumid])
                mean_updraft_w_time[tt,kk] = np.mean(tmp_w[dumid])
                mean_updraft_twc_time[tt,kk] = np.mean(tmp_twc[dumid])

        mean_updraft_tracer_time_avg = np.nanmean(mean_updraft_tracer_time,axis=0)
        mean_updraft_w_time_avg = np.nanmean(mean_updraft_w_time,axis=0)
        mean_updraft_twc_time_avg = np.nanmean(mean_updraft_twc_time,axis=0)
        mean_updraft_tracer_time_avg[mean_updraft_tracer_time_avg == 0.] = np.nan  
        mean_updraft_w_time_avg[mean_updraft_w_time_avg == 0.] = np.nan  
        mean_updraft_twc_time_avg[mean_updraft_twc_time_avg == 0.] = np.nan  
       
    
    
    
    # Average across domain and time
    mean_updraft_tracer = np.zeros(nz)
    mean_updraft_w = np.zeros(nz)
    mean_updraft_twc = np.zeros(nz)
    
    for kk in range(nz):
        dumid = np.where(up_mask_all[:,:,kk,:] == 1.)
        #if (np.size(dumid) == 0.) or (avg_z[kk] < 1000.):
        #    continue
        #else:
        tmp_tracer = tracer_all[:,:,kk,:]
        tmp_w = w_all[:,:,kk,:]
        tmp_twc = twc_all[:,:,kk,:
                         ]
        updraft_tracer = np.ndarray.flatten(tmp_tracer[dumid])
        updraft_w = np.ndarray.flatten(tmp_w[dumid])
        updraft_twc = np.ndarray.flatten(tmp_twc[dumid])
        
        mean_updraft_tracer[kk] = np.nanmean(updraft_tracer)
        mean_updraft_w[kk] = np.nanmean(updraft_w)
        mean_updraft_twc[kk] = np.nanmean(updraft_twc)
        
    mean_updraft_tracer[mean_updraft_tracer == 0.] = np.nan
    mean_updraft_w[mean_updraft_w == 0.] = np.nan
    mean_updraft_twc[mean_updraft_twc == 0.] = np.nan
    
    # Save to diciontary
    #mean_tracer_dict[sim_names[sim_ii]] = mean_updraft_tracer
    #mean_w_dict[sim_names[sim_ii]] = mean_updraft_w
    #mean_twc_dict[sim_names[sim_ii]] = mean_updraft_twc
    #avg_z_dict[sim_names[sim_ii]] = avg_z
    
    
    iplot = False
    if iplot:
        fig = plt.figure(figsize=(16,8))
        Fontsize=20
        ax1 = fig.add_subplot(131)
        ax2 = fig.add_subplot(132)
        ax3 = fig.add_subplot(133)
        axlist = [ax1,ax2,ax3]
        for ax in axlist:
            ax.grid(ls='dotted',lw=1,c='dimgrey')
            ax.tick_params(labelsize=Fontsize)
            ax.set_ylabel('Height [km]',fontsize=Fontsize)
            
        # Tracer
        ax1.set_xlabel('Avg. $\\phi_{updraft}$ [kg$^{-1}$]',fontsize=Fontsize)
        #ax1.plot(mean_updraft_tracer,avg_z*1.e-3,lw=3,c='k')
        ax1.plot(mean_updraft_tracer,z_interps*1.e-3,lw=3,c='k')
        #ax1.plot(mean_updraft_tracer_time_avg,avg_z*1.e-3,lw=3,c='k')
        
        # W
        ax2.set_xlabel('Avg. $w_{updraft}$ [m s$^{-1}$]',fontsize=Fontsize)
        #ax2.plot(mean_updraft_w,avg_z*1.e-3,lw=3,c='k')
        ax2.plot(mean_updraft_w,z_interps*1.e-3,lw=3,c='k')
        #ax2.plot(mean_updraft_w_time_avg,avg_z*1.e-3,lw=3,c='k')

        # TWC
        ax3.set_xlabel('Avg. TWC$_{updraft}$ [g m$^{-3}$]',fontsize=Fontsize)
        #ax3.plot(mean_updraft_twc,avg_z*1.e-3,lw=3,c='k')
        ax3.plot(mean_updraft_twc,z_interps*1.e-3,lw=3,c='k')
        #ax3.plot(mean_updraft_twc_time_avg,avg_z*1.e-3,lw=3,c='k')
        
        plt.show()
        plt.close()
            
    out_dict = {'z':z_interps,\
               'mean_updraft_tracer':mean_updraft_tracer,\
               'mean_updraft_w':mean_updraft_w,\
               'mean_updraft_twc':mean_updraft_twc,\
               }
    
    # Save to dictionary for each simulation
    print('Writing to dictionary...')
    savdir = '/glade/scratch/mckenna/AMIE/amie_post/updraft_tracer/'
    f = open(savdir+sim_names[sim_ii]+'_updraft_tracer.p','wb')
    pickle.dump(out_dict,f)
    f.close()
    print('Completed writing to dictionary. Done.')    
    
    
    
