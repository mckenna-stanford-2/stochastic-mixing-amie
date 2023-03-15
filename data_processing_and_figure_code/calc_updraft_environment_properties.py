#==========================================================================
# Title: calc_updraft_environment_properties.py
# Author: McKenna W. Stanford
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
import pickle
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

sim_names = sim_names[0:2]
num_sims = len(sim_names)

substring = 'stoch'

for sim_ii in range(num_sims):
    print('Simulation:',sim_names[sim_ii])
 
    updraft_qv_all_times = []
    updraft_w_all_times = []
    updraft_tracer_all_times = []
    updraft_mf_all_times = []
    updraft_xkmh_all_times = []
    updraft_xkhh_all_times = []
    updraft_qccon_all_times = []
    updraft_qcevp_all_times = []
    updraft_cconevp_all_times = []
    updraft_mse_all_times = []
    updraft_mse_li_all_times = []
    if substring in sim_names[sim_ii]:
        updraft_rand_pert_all_times = []
    updraft_buoy_all_times = []

    env_qv_all_times = []
    env_tracer_all_times = []
    env_mf_all_times = []
    env_xkmh_all_times = []
    env_xkhh_all_times = []
    env_qccon_all_times = []
    env_qcevp_all_times = []
    env_cconevp_all_times = []
    env_mse_all_times = []
    env_mse_li_all_times = []
    env_buoy_all_times = []


    path = '/glade/scratch/mckenna/AMIE/'+sim_names[sim_ii]+'/'
    files = sorted(glob.glob(path+'wrfout*'))[60:]
    nt = len(files)

    nt = 2
    #========================================================
    # Loop through times
    #========================================================
    for tt in range(nt):
        print('Time step {}/{}'.format(tt+1,nt))
        tmpfile = files[tt]

        
        # Read in with xarray
        ncfile = xarray.open_dataset(tmpfile)
        lat = ncfile['XLAT'].values.squeeze().T
        lon = ncfile['XLONG'].values.squeeze().T
        tracer = ncfile['TR1'].values.squeeze().T
        lon = np.array(ncfile['XLONG'].copy())
        tracer = np.array(ncfile['TR1'].copy())
        xkmh = np.array(ncfile['XKMH'].copy())
        xkhh = np.array(ncfile['XKHH'].copy())
        if substring in sim_names[sim_ii]:
            rand_pert = np.array(ncfile['RAND_PERT'].copy())
        ncfile.close()
        
        # Read in with wrf.getvar
        ncfile = Dataset(tmpfile)
        w = wrf.getvar(ncfile,'wa',meta=False).data.T
        qi = wrf.getvar(ncfile,'QICE',meta=False).data.T
        qr = wrf.getvar(ncfile,'QRAIN',meta=False).data.T
        qc = wrf.getvar(ncfile,'QCLOUD',meta=False).data.T
        qv = wrf.getvar(ncfile,'QVAPOR',meta=False).data.T
        cconevp = wrf.getvar(ncfile,'CCONEVP',meta=False).data.T
        qcevp = wrf.getvar(ncfile,'QCEVP',meta=False).data.T
        qccon = wrf.getvar(ncfile,'QCCON',meta=False).data.T
        tv = wrf.getvar(ncfile,'tv',units='K',meta=False).data.T
        t = wrf.getvar(ncfile,'temp',units='K',meta=False).data.T
        pres = wrf.getvar(ncfile,'pres',units='hPa',meta=False).data.T
        z = wrf.getvar(ncfile,'z',meta=False).data.T
        time = pd.to_datetime(wrf.getvar(ncfile,'times',meta=False))
        ncfile.close()
        
        # Calculate things
        qt = qi+qr+qc
        rd = 287.04
        rho_air = 100.*pres/(rd*tv)
        twc = rho_air*qt*1.e3
        Lv = 2.501*1.e6 #J/kg
        Ls = 2.834*1.e6 #J/kg
        Lf = Ls-Lv
        cp = 1005.7 #J/K/kg
        g = 9.81 #m/s^2
        mse = (t*cp) + (g*z) + (Lv*qv)
        mse_li = (t*cp) + (g*z) + (qv*Lv) - (qi*Lf)
        
        mf = w*rho_air
        
        nx = len(tracer[:,0,0])
        ny = len(tracer[0,:,0])
        nz = len(tracer[0,0,:])
        
        buoy = np.zeros((nx,ny,nz))
        avgz = np.mean(z,axis=(0,1))
        
        avgtv = np.mean(tv,axis=(0,1))
        buoy = np.zeros((nx,ny,nz))
        for kk in range(nz):
            buoy = g*((tv-avgtv)/(avgtv))

        print(aaaaa)

        dx_km = 3.
        dy_km = 3.
        #up_mask = np.zeros((nx,ny,nz))

        updraft_qv_all = []
        updraft_w_all = []
        updraft_tracer_all = []
        updraft_mf_all = []
        updraft_xkmh_all = []
        updraft_xkhh_all = []
        updraft_qccon_all = []
        updraft_qcevp_all = []
        updraft_cconevp_all = []
        updraft_mse_all = []
        updraft_mse_li_all = []
        updraft_buoy_all = []
        if substring in sim_names[sim_ii]:
            updraft_rand_pert_all = []
            
        env_qv_all = []
        env_tracer_all = []
        env_mf_all = []
        env_xkmh_all = []
        env_xkhh_all = []
        env_qccon_all = []
        env_qcevp_all = []
        env_cconevp_all = []
        env_mse_all_times = []
        env_mse_li_all = []
        env_buoy_all = []        
        
        
        
        for kk in range(nz):
            if avgz[kk] < 1000.:
                continue
            avgz_all.append(avgz[kk])
            
            tmp_w = np.array(w[:,:,kk])
            tmp_twc = np.array(twc[:,:,kk])
            tmp_tracer = np.array(tracer[:,:,kk])
            tmp_qv = np.array(qv[:,:,kk])
            tmp_mf = np.array(mf[:,:,kk])
            tmp_xkmh = np.array(xkmh[:,:,kk])
            tmp_xkhh = np.array(xkhh[:,:,kk])
            tmp_qccon = np.array(qccon[:,:,kk])
            tmp_qcevp = np.array(qcevp[:,:,kk])
            tmp_cconevp = np.array(cconevp[:,:,kk])
            tmp_mse = np.array(mse[:,:,kk])
            tmp_mse_li = np.array(mse_li[:,:,kk])
            tmp_buoy = np.array(buoy[:,:,kk])
            if substring in sim_names[sim_ii]:
                tmp_rand_pert = np.array(rand_pert[:,:,kk])
            
            up_mask = np.where((tmp_w > 1.) & (tmp_twc > 0.1))
            num_points = len(up_mask[0,:])
            
            updraft_tracer = []
            env_tracer = []
            updraft_w = []    
            updraft_qv = []
            env_qv = []
            updraft_mf = []
            env_mf = []
            updraft_xkmh = []
            updraft_xkhh = []
            env_xkmh = []
            env_xkhh = []
            updraft_qccon = []
            env_qccon = []
            updraft_qcevp = []
            env_qcevp = []
            updraft_cconevp = []
            env_cconevp = []
            updraft_mse = []
            updraft_mse_li = []
            env_mse = []
            env_mse_li = []
            updraft_buoy = []
            env_buoy = []
            
            if substring in sim_names[sim_ii]:
                updraft_rand_pert = []
            # now loop through number of points on the level
            # identified as an updraft
            
            if num_points > 0.:
                #print('num_points:',num_points)
                for hh in range(num_points):
                    tmpi = up_mask[0,hh]
                    tmpj = up_mask[1,hh]
                    #print(tmpi,tmpj,tmp_w[tmpi,tmpj],tmp_tracer[tmpi,tmpj])
                    updraft_w.append(tmp_w[tmpi,tmpj])
                    updraft_tracer.append(tmp_tracer[tmpi,tmpj])
                    updraft_qv.append(tmp_qv[tmpi,tmpj])
                    updraft_mf.append(tmp_mf[tmpi,tmpj])
                    updraft_xkmh.append(tmp_xkmh[tmpi,tmpj])
                    updraft_xkhh.append(tmp_xkhh[tmpi,tmpj])
                    updraft_qccon.append(tmp_qccon[tmpi,tmpj])
                    updraft_qcevp.append(tmp_qcevp[tmpi,tmpj])
                    updraft_cconevp.append(tmp_cconevp[tmpi,tmpj])
                    updraft_mse.append(tmp_mse[tmpi,tmpj])
                    updraft_mse_li.append(tmp_mse_li[tmpi,tmpj])
                    updraft_buoy.append(tmp_buoy[tmpi,tmpj])
                    if substring in sim:
                        updraft_rand_pert.append(tmp_rand_pert[tmpi,tmpj])
                    #print(aaaa)
                    
                    # now search within 10 km radius of updraft grid point
                    # and calculate environmental tracer
                    bg_thresh = 10
                    bg_thresh_2 = 10
                    i_min = max(tmpi-int(bg_thresh_2/dx_km),0)
                    i_max = min(tmpi+int(bg_thresh_2/dx_km),nx-1)
                    j_min = max(tmpj-int(bg_thresh_2/dx_km),0)
                    j_max = min(tmpj+int(bg_thresh_2/dx_km),ny-1)
                    
                    # create matrix, centered on the updraft, that is size
                    # len_i x len_j (nominally 7 x 7) that will hold the 
                    # radial_distance and the 
                    # tracer values surrounding the updraft. Will
                    # zero out the tracer value on the updraft iteself
                    
                    len_i = i_max-i_min+1
                    len_j = j_max-j_min+1
                    bg_radial_distance = np.zeros((len_i,len_j))
                    bg_tracer = np.zeros((len_i,len_j))
                    bg_qv = np.zeros((len_i,len_j))
                    bg_w = np.zeros((len_i,len_j))
                    bg_mf = np.zeros((len_i,len_j))
                    bg_xkmh = np.zeros((len_i,len_j))
                    bg_xkhh = np.zeros((len_i,len_j))
                    bg_qccon = np.zeros((len_i,len_j))
                    bg_qcevp = np.zeros((len_i,len_j))
                    bg_cconevp = np.zeros((len_i,len_j))
                    bg_mse = np.zeros((len_i,len_j))
                    bg_mse_li = np.zeros((len_i,len_j))
                    bg_buoy = np.zeros((len_i,len_j))
                    
                    dumi = 0
                    for i_bg in range(i_min,(i_max+1)):
                        #print(i_bg)
                        dumj = 0
                        for j_bg in range(j_min,(j_max+1)):
                            #print(i_bg,j_bg)
                            radial_distance = np.sqrt(((dx_km*(tmpi-i_bg))**2.+(dx_km*(tmpj-j_bg))**2.))
                            bg_radial_distance[dumi,dumj] = radial_distance
                            if radial_distance <= bg_thresh:
                                bg_tracer[dumi,dumj] = tmp_tracer[i_bg,j_bg]
                                bg_w[dumi,dumj] = tmp_w[i_bg,j_bg]
                                bg_qv[dumi,dumj] = tmp_qv[i_bg,j_bg]
                                bg_mf[dumi,dumj] = tmp_mf[i_bg,j_bg]
                                bg_xkmh[dumi,dumj] = tmp_xkmh[i_bg,j_bg]
                                bg_xkhh[dumi,dumj] = tmp_xkhh[i_bg,j_bg]
                                bg_qccon[dumi,dumj] = tmp_qccon[i_bg,j_bg]
                                bg_qcevp[dumi,dumj] = tmp_qcevp[i_bg,j_bg]
                                bg_cconevp[dumi,dumj] = tmp_cconevp[i_bg,j_bg]
                                bg_mse[dumi,dumj] = tmp_mse[i_bg,j_bg]
                                bg_mse_li[dumi,dumj] = tmp_mse_li[i_bg,j_bg]
                                bg_buoy[dumi,dumj] = tmp_buoy[i_bg,j_bg]
                            elif radial_distance > bg_thresh:
                                bg_tracer[dumi,dumj] = np.nan
                                bg_qv[dumi,dumj] = np.nan
                                bg_mf[dumi,dumj] = np.nan
                                bg_xkmh[dumi,dumj] = np.nan
                                bg_xkhh[dumi,dumj] = np.nan
                                bg_qccon[dumi,dumj] = np.nan
                                bg_qcevp[dumi,dumj] = np.nan
                                bg_cconevp[dumi,dumj] = np.nan
                                bg_mse[dumi,dumj] = np.nan
                                bg_mse_li[dumi,dumj] = np.nan
                                bg_buoy[dumi,dumj] = np.nan
                                #bg_w[dumi,dumj] = np.nan
                            if radial_distance == 0.:
                                bg_tracer[dumi,dumj] = np.nan
                                bg_qv[dumi,dumj] = np.nan
                                bg_mf[dumi,dumj] = np.nan
                                bg_xkmh[dumi,dumj] = np.nan
                                bg_xkhh[dumi,dumj] = np.nan
                                bg_qccon[dumi,dumj] = np.nan
                                bg_qcevp[dumi,dumj] = np.nan
                                bg_cconevp[dumi,dumj] = np.nan
                                bg_mse[dumi,dumj] = np.nan
                                bg_mse_li[dumi,dumj] = np.nan
                                bg_buoy[dumi,dumj] = np.nan

                            dumj+=1
                        dumi+=1
                        
                    env_tracer.append(np.nanmean(bg_tracer))
                    env_qv.append(np.nanmean(bg_qv))
                    env_mf.append(np.nanmean(bg_mf))
                    env_xkmh.append(np.nanmean(bg_xkmh))
                    env_xkhh.append(np.nanmean(bg_xkhh))
                    env_qccon.append(np.nanmean(bg_qccon))
                    env_qcevp.append(np.nanmean(bg_qcevp))
                    env_cconevp.append(np.nanmean(bg_cconevp))
                    env_mse.append(np.nanmean(bg_mse))
                    env_mse_li.append(np.nanmean(bg_mse_li))
                    env_buoy.append(np.nanmean(bg_buoy))
                    
                    if False:
                        # test plot that countours in lines the radial distance
                        # again in lines the vertical velocity, and in color
                        # the tracer concentration
                        class nf(float):
                            def __repr__(self):
                                s = f'{self:.1f}'
                                return f'{self:.0f}' if s[-1] == '0' else s                    


                        fig = plt.figure(figsize=(10,8))
                        ax = fig.add_subplot(111)
                        CS = ax.contour(bg_radial_distance,colors='black')

                        CS.levels = [nf(val) for val in CS.levels]
                        # Label levels with specially formatted floats
                        if plt.rcParams["text.usetex"]:
                            fmt = r'%r km'
                        else:
                            fmt = '%r km'
                        ax.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=15)
                        ax.set_xlabel('x grid points',fontsize=20)
                        ax.set_ylabel('y grid points',fontsize=20)
                        plt.tick_params(labelsize=20)
                        tmpplot = ax.contourf(bg_w,cmap='seismic',\
                                              levels=np.arange(-3,3.1,0.1))
                        cbar=plt.colorbar(tmpplot,ticks=[-3,-2.5,-2,-1.5,-1,-0.5,0,\
                                                        0,0.5,1,1.5,2,2.5,3])
                        cbar.ax.set_ylabel('w [m s$^{-1}$]',fontsize=20)
                        cbar.ax.tick_params(labelsize=15)

                        tracer_plot = ax.contour(bg_tracer,\
                                    levels=[0.,0.6,0.7,0.8,0.9,1.],\
                                    colors='grey',\
                                    linewidths=4
                                    )

                        tracer_plot.levels = [nf(val) for val in tracer_plot.levels]
                        # Label levels with specially formatted floats
                        if plt.rcParams["text.usetex"]:
                            fmt = r'%r g kg$^{-1}$'
                        else:
                            fmt = '%r g kg$^{-1}$'
                        tmp=ax.clabel(tracer_plot, tracer_plot.levels,\
                                  inline=True, fmt=fmt, fontsize=15)
                        [txt.set_bbox(dict(facecolor='white', \
                            edgecolor='none', pad=0)) for txt in tmp]
                        plt.show()
                        if hh > 12:
                            print(aaaaa)
                            
            updraft_tracer_all.append(updraft_tracer)
            env_tracer_all.append(env_tracer)
            updraft_w_all.append(updraft_w)
            updraft_qv_all.append(updraft_qv)
            env_qv_all.append(env_qv)
            env_mf_all.append(env_mf)
            updraft_mf_all.append(updraft_mf)
            updraft_xkmh_all.append(updraft_xkmh)
            updraft_xkhh_all.append(updraft_xkhh)
            env_xkmh_all.append(env_xkmh)
            env_xkhh_all.append(env_xkhh)
            updraft_qccon_all.append(updraft_qccon)
            updraft_qcevp_all.append(updraft_qcevp)
            updraft_cconevp_all.append(updraft_cconevp)
            env_qccon_all.append(env_qccon)
            env_qcevp_all.append(env_qcevp)
            env_cconevp_all.append(env_cconevp)
            updraft_mse_all.append(updraft_mse)
            updraft_mse_li_all.append(updraft_mse_li)
            env_mse_all.append(env_mse)
            env_mse_li_all.append(env_mse_li)
            env_buoy_all.append(env_buoy)
            if substring in sim:
                updraft_rand_pert_all.append(updraft_rand_pert)
            
        updraft_tracer_all = np.array(updraft_tracer_all)
        updraft_w_all = np.array(updraft_w_all)
        updraft_qv_all = np.array(updraft_qv_all)
        env_tracer_all = np.array(env_tracer_all)
        env_qv_all = np.array(env_qv_all)
        updraft_mf_all = np.array(updraft_mf_all)
        env_mf_all = np.array(env_mf_all)
        updraft_xkmh_all = np.array(updraft_xkmh_all)
        updraft_xkhh_all = np.array(updraft_xkhh_all)
        env_xkmh_all = np.array(env_xkmh_all)
        env_xkhh_all = np.array(env_xkhh_all)
        updraft_qccon_all = np.array(updraft_qccon_all)
        updraft_qcevp_all = np.array(updraft_qcevp_all)
        updraft_cconevp_all = np.array(updraft_cconevp_all)
        env_qccon_all = np.array(env_qccon_all)
        env_qcevp_all = np.array(env_qcevp_all)
        env_cconevp_all = np.array(env_cconevp_all)
        updraft_mse_all = np.array(updraft_mse_all)
        updraft_mse_li_all = np.array(updraft_mse_li_all)
        env_mse_all = np.array(env_mse_all)
        env_mse_li_all = np.array(env_mse_li_all)
        env_buoy_all = np.array(env_buoy_all)
        updraft_buoy_all = np.array(updraft_buoy_all)
        if substring in sim:
            updraft_rand_pert_all = np.array(updraft_rand_pert_all)

        updraft_tracer_all_times.append(updraft_tracer_all)
        env_tracer_all_times.append(env_tracer_all)
        updraft_w_all_times.append(updraft_w_all)
        updraft_qv_all_times.append(updraft_qv_all)
        env_qv_all_times.append(env_qv_all)
        env_mf_all_times.append(env_mf_all)
        updraft_mf_all_times.append(updraft_mf_all)
        updraft_xkmh_all_times.append(updraft_xkmh_all)
        updraft_xkhh_all_times.append(updraft_xkhh_all)
        env_xkmh_all_times.append(env_xkmh_all)
        env_xkhh_all_times.append(env_xkhh_all)
        updraft_qccon_all_times.append(updraft_qccon_all)
        updraft_qcevp_all_times.append(updraft_qcevp_all)
        updraft_cconevp_all_times.append(updraft_cconevp_all)
        env_qccon_all_times.append(env_qccon_all)
        env_qcevp_all_times.append(env_qcevp_all)
        env_cconevp_all_times.append(env_cconevp_all)
        updraft_mse_all_times.append(updraft_mse_all)
        updraft_mse_li_all_times.append(updraft_mse_li_all)
        env_mse_all_times.append(env_mse_all)
        env_mse_li_all_times.append(env_mse_li_all)
        env_buoy_all_times.append(env_buoy_all)
        updraft_buoy_all_times.append(updraft_buoy_all)
        if substring in sim:
            updraft_rand_pert_all_times.append(updraft_rand_pert_all)            

            
    updraft_qv_all_times = np.array(updraft_qv_all_times)
    updraft_w_all_times = np.array(updraft_w_all_times)
    updraft_tracer_all_times = np.array(updraft_tracer_all_times)
    env_tracer_all_times = np.array(env_tracer_all_times)
    env_qv_all_times = np.array(env_qv_all_times)
    updraft_mf_all_times = np.array(updraft_mf_all_times)
    env_mf_all_times = np.array(env_mf_all_times)
    updraft_xkmh_all_times = np.array(updraft_xkmh_all_times)
    updraft_xkhh_all_times = np.array(updraft_xkhh_all_times)
    env_xkmh_all_times = np.array(env_xkmh_all_times)
    env_xkhh_all_times = np.array(env_xkhh_all_times)
    updraft_qccon_all_times = np.array(updraft_qccon_all_times)
    updraft_qcevp_all_times = np.array(updraft_qcevp_all_times)
    updraft_cconevp_all_times = np.array(updraft_cconevp_all_times)
    updraft_mse_all_times = np.array(updraft_mse_all_times)
    updraft_mse_li_all_times = np.array(updraft_mse_li_all_times)
    env_qccon_all_times = np.array(env_qccon_all_times)
    env_qcevp_all_times = np.array(env_qcevp_all_times)
    env_cconevp_all_times = np.array(env_cconevp_all_times)
    env_mse_all_times = np.array(env_mse_all_times)
    env_mse_li_all_times = np.array(env_mse_li_all_times)
    if substring in sim:
        updraft_rand_pert_all_times = np.array(updraft_rand_pert_all_times)
    updraft_buoy_all_times = np.array(updraft_buoy_all_times)
    env_buoy_all_times = np.array(env_buoy_all_times)
            
     
    if substring in sim:
        up_env_dict = {'updraft_tracer':updraft_tracer_all_times,\
                         'updraft_w':updraft_w_all_times,\
                         'updraft_qv':updraft_qv_all_times,\
                         'env_tracer':env_tracer_all_times,\
                         'env_qv':env_qv_all_times,\
                         'avgz':avgz,\
                         'time':time,\
                         'updraft_mf':updraft_mf_all_times,\
                         'env_mf':env_mf_all_times,\
                         'updraft_xkmh':updraft_xkmh_all_times,\
                         'updraft_xkhh':updraft_xkhh_all_times,\
                         'env_xkmh':env_xkmh_all_times,\
                         'env_xkhh':env_xkhh_all_times,\
                         'updraft_qccon':updraft_qccon_all_times,\
                         'updraft_qcevp':updraft_qcevp_all_times,\
                         'updraft_cconevp':updraft_cconevp_all_times,\
                         'env_qccon':env_qccon_all_times,\
                         'env_qcevp':env_qcevp_all_times,\
                         'env_cconevp':env_cconevp_all_times,\
                         'updraft_mse':updraft_mse_all_times,\
                         'updraft_mse_li':updraft_mse_li_all_times,\
                         'env_mse':env_mse_all_times,\
                         'env_mse_li':env_mse_li_all_times,\
                         'updraft_rand_pert':updraft_rand_pert_all_times,\
                         'updraft_buoy':updraft_buoy_all_times,\
                         'env_buoy':env_buoy_all_times,\
                        }
    else:
        up_env_dict = {'updraft_tracer':updraft_tracer_all_times,\
                         'updraft_w':updraft_w_all_times,\
                         'updraft_qv':updraft_qv_all_times,\
                         'env_tracer':env_tracer_all_times,\
                         'env_qv':env_qv_all_times,\
                         'avgz':avgz,\
                         'time':time,\
                         'updraft_mf':updraft_mf_all_times,\
                         'env_mf':env_mf_all_times,\
                         'updraft_xkmh':updraft_xkmh_all_times,\
                         'updraft_xkhh':updraft_xkhh_all_times,\
                         'env_xkmh':env_xkmh_all_times,\
                         'env_xkhh':env_xkhh_all_times,\
                         'updraft_qccon':updraft_qccon_all_times,\
                         'updraft_qcevp':updraft_qcevp_all_times,\
                         'updraft_cconevp':updraft_cconevp_all_times,\
                         'env_qccon':env_qccon_all_times,\
                         'env_qcevp':env_qcevp_all_times,\
                         'env_cconevp':env_cconevp_all_times,\
                         'updraft_mse':updraft_mse_all_times,\
                         'updraft_mse_li':updraft_mse_li_all_times,\
                         'env_mse':env_mse_all_times,\
                         'env_mse_li':env_mse_li_all_times,\
                         'updraft_buoy':updraft_buoy_all_times,\
                         'env_buoy':env_buoy_all_times,\
                        }                    
        
        
        
    savdir = '/glade/scratch/mckenna/AMIE/'+sim_names[sim_ii]+'/'
    f = open(savdir+sim_names[sim_ii]+'_updraft_environment_dict.p','wb')
    pickle.dump(up_env_dict,f)
    f.close()
        
    
#print(aaaaa)

print('done')

