#==========================================================================
# Title: coarse_grain_spol.py
# Author: McKenna W. Stanford
# Utility: Coarse-grains 1-km simulation to 3-km grid spacing.
# Writes this and native res. variables to a dictionary named 
# spol_coarse_grained.p.
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
import pickle
import pandas as pd
import warnings
warnings.filterwarnings("ignore", message="RuntimeWarning: Mean of empty slice")

# read in spol reflectivity
path = '/glade/scratch/mckenna/AMIE/spol/refl/'
refl_files = sorted(glob.glob(path+'*.nc'))[44:]
nt = len(refl_files)

# read in spol rain rates
path = '/glade/scratch/mckenna/AMIE/spol/rain_rates/'
rain_rate_files = sorted(glob.glob(path+'*.nc'))[44:]


dbz_deg_arr = []
rain_rate_deg_arr = []
rain_rate_min_deg_arr = []
rain_rate_max_deg_arr = []
time_arr = []

dbz_arr = []
rain_rate_arr = []
rain_rate_min_arr = []
rain_rate_max_arr = []

for tt in range(nt):
    print('Time step {}/{}:'.format(tt+1,nt))

    # Read in reflectivity file
    ncfile = xarray.open_dataset(refl_files[tt])
    lat = ncfile['lat0'].values.squeeze().T
    lon = ncfile['lon0'].values.squeeze().T
    dbz = ncfile['REFL'].values.squeeze().T
    z = ncfile['z0'].values.squeeze() # km
    time = pd.to_datetime(ncfile['time'].values[0])
    ncfile.close()
    
    # Read in rain rates
    ncfile = xarray.open_dataset(rain_rate_files[tt])
    rain_rate = ncfile['rain_rate'].values.squeeze().T
    rain_rate_max = ncfile['rain_rate_max'].values.squeeze().T
    rain_rate_min = ncfile['rain_rate_min'].values.squeeze().T
    ncfile.close()


    nz = len(z)
    nx = 100
    ny = 100
    dbz_cg = np.zeros((nx,ny,nz))-999.
    coarseness = 3
    
    # Coase grain rain rates
    tmp_lon = lon.reshape((lon.shape[0] // coarseness, coarseness,
                            lon.shape[1] // coarseness, coarseness))
    coarse_lon = np.nanmean(tmp_lon, axis=(1,3))   
    
    tmp_lat = lat.reshape((lat.shape[0] // coarseness, coarseness,
                            lat.shape[1] // coarseness, coarseness))
    coarse_lat = np.nanmean(tmp_lat, axis=(1,3))           
        
    tmp_rain_rate = rain_rate.reshape((rain_rate.shape[0] // coarseness, coarseness,\
                             rain_rate.shape[1] // coarseness, coarseness))
    tmp_rain_rate_min = rain_rate_min.reshape((rain_rate_min.shape[0] // coarseness, coarseness,\
                             rain_rate_min.shape[1] // coarseness, coarseness))
    tmp_rain_rate_max = rain_rate_max.reshape((rain_rate_max.shape[0] // coarseness, coarseness,\
                             rain_rate_max.shape[1] // coarseness, coarseness)) 
    
    coarse_rain_rate = np.nanmean(tmp_rain_rate,axis=(1,3))
    coarse_rain_rate_min = np.nanmean(tmp_rain_rate_min,axis=(1,3))
    coarse_rain_rate_max = np.nanmean(tmp_rain_rate_max,axis=(1,3))
    
    iplot = False
    if iplot:
        rain_rate_levs = 10.**np.arange(-2,2.1,0.1)
        fig = plt.figure(figsize=(10,4))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        dum_ticks = [1.e-2,1.e-1,1.e0,1.e1,1.e2]
        plot1=ax1.contourf(lon,lat,rain_rate,levels=rain_rate_levs,\
                           cmap='cividis',norm=matplotlib.colors.LogNorm(),extend='both')
        cbar1 = fig.colorbar(plot1,ax=ax1,ticks=dum_ticks)
        cbar1.ax.set_ylabel('Rain Rate [mm hr$^{-1}$]')
        
        plot2=ax2.contourf(coarse_lon,coarse_lat,coarse_rain_rate,levels=rain_rate_levs,\
                           cmap='cividis',norm=matplotlib.colors.LogNorm(),extend='both')
        cbar2 = fig.colorbar(plot2,ax=ax2,ticks=dum_ticks)
        cbar2.ax.set_ylabel('Rain Rate [mm hr$^{-1}$]')
        
        ax1.set_title('Native')
        ax2.set_title('Coarse-Grained')
        plt.show()
        plt.close()
            
        
        
    # Coarse-grain reflectivity at all heights
    for kk in range(0,nz):
        #print('Height:',z[kk])
        tmp_dbz = dbz[:,:,kk]
        ze = 10.**(tmp_dbz/10.)        
        tmp_ze = ze.reshape((ze.shape[0] // coarseness, coarseness,
                            ze.shape[1] // coarseness, coarseness))
        
        coarse_ze = np.nanmean(tmp_ze, axis=(1,3))  
        coarse_dbz = 10.*np.log10(coarse_ze)
        dbz_cg[:,:,kk] = coarse_dbz

        
        
    iplot = False
    if iplot:
        dbz_levs = np.arange(0,61,1)
        fig = plt.figure(figsize=(10,4))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        dum_ticks = np.arange(0,70,10)
        plot1=ax1.contourf(lon,lat,dbz[:,:,4],levels=dbz_levs,\
                           cmap='cividis')
        cbar1 = fig.colorbar(plot1,ax=ax1,ticks=dum_ticks)
        cbar1.ax.set_ylabel('Refl [dBZ]')
        
        plot2=ax2.contourf(coarse_lon,coarse_lat,dbz_cg[:,:,4],levels=dbz_levs,\
                           cmap='cividis')
        cbar2 = fig.colorbar(plot2,ax=ax2,ticks=dum_ticks)
        cbar2.ax.set_ylabel('Refl [dBZ]')
        
        ax1.set_title('Native')
        ax2.set_title('Coarse-Grained')
        plt.show()
        plt.close()
        
        
    
    time_arr.append(time)    
    rain_rate_arr.append(rain_rate)
    rain_rate_max_arr.append(rain_rate_max)
    rain_rate_min_arr.append(rain_rate_min)
    dbz_arr.append(dbz)
    dbz_deg_arr.append(dbz_cg)
    rain_rate_deg_arr.append(coarse_rain_rate)
    rain_rate_max_deg_arr.append(coarse_rain_rate_max)
    rain_rate_min_deg_arr.append(coarse_rain_rate_min)
    
time = np.array(time)
rain_rate_arr = np.array(rain_rate_arr)
rain_rate_max_arr = np.array(rain_rate_max_arr)
rain_rate_min_arr = np.array(rain_rate_min_arr)
dbz_arr = np.array(dbz_arr)
rain_rate_deg_arr = np.array(rain_rate_deg_arr)
rain_rate_max_deg_arr = np.array(rain_rate_max_deg_arr)
rain_rate_min_deg_arr = np.array(rain_rate_min_deg_arr)
dbz_deg_arr = np.array(dbz_deg_arr)


# Write to dictionary
spol_dict = {'1km':{'lon':lon,\
                    'lat':lat,\
                    'z':z,\
                    'time':time_arr,\
                    'rain_rate':rain_rate_arr,\
                    'rain_rate_min':rain_rate_min_arr,\
                    'rain_rate_max':rain_rate_max_arr,\
                    'dbz':dbz_arr,\
                    },\
             '3km':{'lon':coarse_lon,\
                    'lat':coarse_lat,\
                    'z':z,\
                    'time':time_arr,\
                    'rain_rate':rain_rate_deg_arr,\
                    'rain_rate_min':rain_rate_min_deg_arr,\
                    'rain_rate_max':rain_rate_max_deg_arr,\
                    'dbz':dbz_deg_arr,\
                    }
            }

# Print diagnostics
#if False:
if True:
    for key,val in spol_dict.items():
        print(key)
        
        for key2,val2 in spol_dict[key].items():
            print('   ',key2,np.shape(val2))
    
print('Writing to dictionary...')
savdir = '/glade/scratch/mckenna/AMIE/amie_post/spol/'
f = open(savdir+'spol_coarse_grained.p','wb')
pickle.dump(spol_dict,f)
f.close()
        
print('Done.')
        

        
