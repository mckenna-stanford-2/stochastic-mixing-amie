#==========================================================================
# Title: compute_wrf_int_mass_flux_1km_CG.py
# Author: McKenna W. Stanford
# Utility: Vertically integrates mass flux considering only upward or downward motions
# and upward or downward motions conditionally sampled for grid points with
# vertical velocity > 1 m/s.
# Writes to a dictionary named '<SIM_NAME>_int_mass_flux_dict.p'
#==========================================================================

#------------------------------------------
# Imports
#------------------------------------------
import numpy as np
import glob
import xarray
import matplotlib
import pickle
import pandas as pd
from netCDF4 import Dataset
import wrf

path = '/glade/scratch/mckenna/AMIE/amie_post/mass_flux/'

file = path+'baseline_1km_CG_mass_flux_twc_w_dict.p'
pkl_file = open(file,'rb')
pkl_dict = pickle.load(pkl_file)
nt = len(pkl_dict['time'])
nz = len(pkl_dict['z'][0,0,0,:])
nx = len(pkl_dict['z'][0,:,0,0])
ny = len(pkl_dict['z'][0,0,:,0])

# Vertically mass flux and convective mass flux
print('Vertically integating mass flux...')
print('Looping through times...')

mf_int_pos_up = np.zeros((nx,ny,nt))
mf_int_neg_down = np.zeros((nx,ny,nt))
mf_int_conv_up = np.zeros((nx,ny,nt))

for tt in range(nt):
    print('Time step {}/{}'.format(tt+1,nt))
    tmp_z = pkl_dict['z'][tt,:,:,:]
    tmp_w = pkl_dict['w'][tt,:,:,:]
    tmp_mf = pkl_dict['mass_flux'][tt,:,:,:]
    tmp_twc = pkl_dict['twc'][tt,:,:,:]
    nz = len(tmp_z[0,0,:])
    nx = len(tmp_z[:,0,0])
    ny = len(tmp_z[0,:,0])

    #----------------------------------------------
    # Positive upward vertically-integrated
    # mass flux
    #----------------------------------------------
    pos_up_mf = np.zeros((nx,ny,nz))
    pos_up_id = np.where(tmp_mf > 0.)
    pos_up_mf[pos_up_id] = tmp_mf[pos_up_id]
    pos_up_int_mf = np.trapz(pos_up_mf,tmp_z,axis=2)
    mf_int_pos_up[:,:,tt] = pos_up_int_mf     

    #----------------------------------------------
    # Negative downward vertically-integrated
    # mass flux
    #----------------------------------------------
    neg_down_mf = np.zeros((nx,ny,nz))
    neg_down_id = np.where(tmp_mf < 0.)
    neg_down_mf[neg_down_id] = tmp_mf[neg_down_id]
    neg_down_int_mf = np.trapz(neg_down_mf,tmp_z,axis=2)
    mf_int_neg_down[:,:,tt] = neg_down_int_mf     

    #----------------------------------------------
    # Convective updraft mass flux
    #----------------------------------------------
    conv_up_mf = np.zeros((nx,ny,nz))
    conv_up_id = np.where( (tmp_w > 1.) & (tmp_twc > 0.1) )
    conv_up_mf[conv_up_id] = tmp_mf[conv_up_id]
    conv_up_int_mf = np.trapz(conv_up_mf,tmp_z,axis=2)
    mf_int_conv_up[:,:,tt] = conv_up_int_mf     


out_dict = {'pos_up_mf':mf_int_pos_up,\
            'neg_down_mf':mf_int_neg_down,\
            'conv_up_mf':mf_int_conv_up,\
            'time':pkl_dict['time'],\
            'lat':pkl_dict['lat'],\
            'lon':pkl_dict['lon'],\
           }

savdir = '/glade/scratch/mckenna/AMIE/amie_post/mass_flux/'
f = open(savdir+'baseline_1km_CG_int_mass_flux_dict.p','wb')
pickle.dump(out_dict,f)
f.close()
print('Completed writing to dictionary. Done.')
