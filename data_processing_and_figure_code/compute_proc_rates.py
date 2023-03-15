#==========================================================================
# Title: compute_proc_rates.py
# Author: McKenna W. Stanford
# Utility: Reads in wrfout* files and grabs various process rates and subtracts
# them from the previous file in time, such that the output is accumulated 
# process rates per 15-minutes (outfile time step). Then computes total
# evaporation and condensation by summing across hydrometeor species for
# source/sink processes. Then vertically integrates the total evaporation/
# condensation.
# Writes this to a dictionary named '<SIM_NAME>_proc_rates.p.
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
import pandas as pd
from matplotlib.lines import Line2D

#plt.rc('text',usetex=True)

fig_path = '/glade/u/home/mckenna/figures/'

path = '/glade/scratch/mckenna/AMIE/'
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
            ]
#sim_names = sim_names[0:8]
sim_names = sim_names[8:]
num_sims = len(sim_names)


# Loop through simulations
for sim_ii in range(num_sims):
    print('Simulation:',sim_names[sim_ii])
    files = sorted(glob.glob(path+sim_names[sim_ii]+'/wrfout*'))
    files = files[47:]
    nt = len(files) # number of files/times
    #nt = 3
    
    # Loop through times
    rainnc_arr = []
    snownc_arr = []
    qrevp_arr = []
    qrcon_arr = []
    qcevp_arr = []
    qccon_arr = []
    qidep_arr = []
    qisub_arr = []
    z_arr = []
    time_arr = []
    rho_air_arr = []
    
    for tt in range(nt):
        print('Time step {}/{}:'.format(tt+1,nt))
        
        # Use xarray
        ncfile = xarray.open_dataset(files[tt])
        rainnc = ncfile['RAINNC'].values.squeeze().T
        snownc = ncfile['SNOWNC'].values.squeeze().T
        qrevp = ncfile['QREVP'].values.squeeze().T
        qrcon = ncfile['QRCON'].values.squeeze().T
        qcevp = ncfile['QCEVP'].values.squeeze().T
        qccon = ncfile['QCCON'].values.squeeze().T
        qidep = ncfile['QIDEP'].values.squeeze().T
        qisub = ncfile['QISUB'].values.squeeze().T
        if tt == 0.:
            lat = ncfile['XLAT'].values.squeeze().T
            lon = ncfile['XLONG'].values.squeeze().T  
            lon_1d = lon[:,0]
            lat_1d = lat[0,:]
        time = pd.to_datetime(ncfile['XTIME'].values[0])
        ncfile.close()
        
        time_arr.append(time)
        rainnc_arr.append(rainnc)
        snownc_arr.append(snownc)
        qrevp_arr.append(qrevp)
        qrcon_arr.append(qrcon)
        qcevp_arr.append(qcevp)
        qccon_arr.append(qccon)
        qidep_arr.append(qidep)
        qisub_arr.append(qisub)
        
        # Use wrf python package
        ncfile = Dataset(files[tt])
        z = wrf.getvar(ncfile,'z',meta=False).data.T # height AGL, units of m
        tv = wrf.getvar(ncfile,"tv",units='K',meta=False).data.T # virtual temperature in K
        pressure = wrf.getvar(ncfile,"pres",units='hPa',meta=False).data.T # pressure in hPa
        ncfile.close()
        
        rd = 287.04
        rho_air = 100.*pressure/(rd*tv) #kg/m^3 

        
        z_arr.append(z)
        rho_air_arr.append(rho_air)

    rainnc_arr = np.array(rainnc_arr)
    snownc_arr = np.array(snownc_arr)
    qrevp_arr = np.array(qrevp_arr)
    qrcon_arr = np.array(qrcon_arr)
    qcevp_arr = np.array(qcevp_arr)
    qccon_arr = np.array(qccon_arr)
    qidep_arr = np.array(qidep_arr)
    qisub_arr = np.array(qisub_arr)
    z_arr = np.array(z_arr)
    rho_air_arr = np.array(rho_air_arr)
    time_arr = np.array(time_arr)
    
    # RAINNC is accumualted in time for each output file
    # Need to subtract to get values per time dimension
    # So, these will have values of mm/(15 minutes)
    rainnc_arr_sub = []
    snownc_arr_sub = []
    qrevp_arr_sub = []
    qrcon_arr_sub = []
    qcevp_arr_sub = []
    qccon_arr_sub = []
    qidep_arr_sub = []
    qisub_arr_sub = []
    
    print('Subtracting along time axis...')
    for tt in range(nt-1):
        print('time step {}/{}:'.format(tt+1,nt-1))

        rainnc_arr_sub.append(rainnc_arr[tt+1,:,:]-rainnc_arr[tt,:,:])
        snownc_arr_sub.append(snownc_arr[tt+1,:,:]-snownc_arr[tt,:,:])
        qrevp_arr_sub.append(qrevp_arr[tt+1,:,:]-qrevp_arr[tt,:,:])
        qrcon_arr_sub.append(qrcon_arr[tt+1,:,:]-qrcon_arr[tt,:,:])
        qcevp_arr_sub.append(qcevp_arr[tt+1,:,:]-qcevp_arr[tt,:,:])
        qccon_arr_sub.append(qccon_arr[tt+1,:,:]-qccon_arr[tt,:,:])
        qidep_arr_sub.append(qidep_arr[tt+1,:,:]-qidep_arr[tt,:,:])
        qisub_arr_sub.append(qisub_arr[tt+1,:,:]-qisub_arr[tt,:,:])
        
    rainnc_arr_sub = np.array(rainnc_arr_sub)
    snownc_arr_sub = np.array(snownc_arr_sub)
    qrevp_arr_sub = np.array(qrevp_arr_sub)
    qrcon_arr_sub = np.array(qrcon_arr_sub)
    qcevp_arr_sub = np.array(qcevp_arr_sub)
    qccon_arr_sub = np.array(qccon_arr_sub)
    qidep_arr_sub = np.array(qidep_arr_sub)
    qisub_arr_sub = np.array(qisub_arr_sub)
    
    # Limit time,z, and rho_air arrays to the time after the first time
    # step to match with subtractions
    time_arr = time_arr[1:]
    z_arr = z_arr[1:,:,:,:]
    rho_air_arr = rho_air_arr[1:,:,:,:]

    # Compute total condensation and evaporation
    print('Computing total condensation and evaporation...')
    tot_cond = qrcon_arr_sub + qccon_arr_sub + qidep_arr_sub
    tot_evap = qrevp_arr_sub + qcevp_arr_sub + qisub_arr_sub
    print('Completed computing total condensation and evaporation.')

    
    # Vertically integrate
    print('Vertically integating total condensation and evaporation vertically...')
    col_int_cond = np.trapz(tot_cond*rho_air_arr,z_arr,axis=3)
    col_int_evap = np.trapz(tot_evap*rho_air_arr,z_arr,axis=3)
    print('Completed vertically integating total condensation and evaporation vertically.')
    
    print('Writing to dictionary...')
    # Put in dictionary and write to pickle file
    out_dict = {'rainnc':rainnc_arr_sub,\
                'snownc':snownc_arr_sub,\
                'qrevp':qrevp,\
                'qrcon':qrcon_arr_sub,\
                'qcevp':qcevp_arr_sub,\
                'qccon':qccon,\
                'qidep':qidep_arr_sub,\
                'qisub':qisub_arr_sub,\
                'rho_air':qisub_arr_sub,\
                'tot_cond':tot_cond,\
                'tot_evap':tot_evap,\
                'col_int_cond':col_int_cond,\
                'col_int_evap':col_int_evap,\
                'z':z_arr,\
                'time':time_arr,\
                'lat':lat,\
                'lon':lon,\
                'lat_1d':lat_1d,\
                'lon_1d':lon_1d}
        
    # Print diagnostics
    if False:
    #if True:
        for key,val in out_dict.items():
            if key != 'time':
                print(key,np.shape(val))
            else:
                print(key,val)
            
            
    savdir = '/glade/scratch/mckenna/AMIE/amie_post/proc_rates/'
    f = open(savdir+sim_names[sim_ii]+'_proc_rates.p','wb')
    pickle.dump(out_dict,f)
    f.close()
    print('Completed writing to dictionary. Done.')


