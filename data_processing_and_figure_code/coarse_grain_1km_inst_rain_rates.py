#==========================================================================
# Title: coarse_grain_1km_inst_rain_rates.py
# Author: McKenna Stanford
# Utility: Coarse-grains 1-km rain rates to 3-km grid spacing.
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
import pickle
from netCDF4 import Dataset
import wrf
from scipy.io import netcdf

path = '/glade/scratch/mckenna/AMIE/amie_post/inst_rain_rates/'

pkl_file = open(path+'baseline_1km_inst_rain_rates.p','rb')
tmpfile = pickle.load(pkl_file)
pkl_file.close()

rain_rate = tmpfile['rain_rate']
lon = tmpfile['lon']
lat = tmpfile['lat']
time = tmpfile['time']
nt = len(time)
nx = len(rain_rate[0,:,0])
ny = len(rain_rate[0,0,:])

rain_rate_3km_all_times = []

for tt in range(nt):
    print('Time step {}/{}:'.format(tt+1,nt))  

    fac = 3.
    nx_3km = int(np.floor(nx/fac))
    ny_3km = int(np.floor(ny/fac))
    

    rain_rate_3km = np.zeros((nx_3km, ny_3km))
    lon_3km = np.zeros((nx_3km, ny_3km))
    lat_3km = np.zeros((nx_3km, ny_3km))
    
    fac = int(fac)
    for ii in range(nx_3km):
        for jj in range(ny_3km):
            rain_rate_3km[ii,jj] = np.mean(rain_rate[tt,ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac)])
            lon_3km[ii,jj] = np.mean(lon[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac)])
            lat_3km[ii,jj] = np.mean(lat[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac)])
    rain_rate_3km_all_times.append(rain_rate_3km)
rain_rate_3km_all_times = np.array(rain_rate_3km_all_times)
             
#levs = [0.01,0.05,0.1,0.5,1,5,10,50,100]
#plt.contourf(lon_3km,lat_3km,rr_3km_all_times[-1,:,:,1],levels=levs,cmap='inferno')

    
rain_rate_deg_dict = {'rain_rate':rain_rate_3km_all_times,\
                      'lon':lon_3km,\
                      'lat':lat_3km,\
                      'time':time}

for key,val in rain_rate_deg_dict.items():
    print(key,np.shape(val))
    
    
iplot = False
if iplot:
    rain_rate_levs = 10.**np.arange(-2,2.1,0.1)
    fig = plt.figure(figsize=(10,4))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    dum_ticks = [1.e-2,1.e-1,1.e0,1.e1,1.e2]
    plot1=ax1.contourf(lon,lat,rain_rate[-1,:,:],levels=rain_rate_levs,\
                       cmap='cividis',norm=matplotlib.colors.LogNorm(),extend='both')
    cbar1 = fig.colorbar(plot1,ax=ax1,ticks=dum_ticks)
    cbar1.ax.set_ylabel('Rain Rate [mm hr$^{-1}$]')

    plot2=ax2.contourf(lon_3km,lat_3km,rain_rate_deg_dict['rain_rate'][-1,:,:],levels=rain_rate_levs,\
                       cmap='cividis',norm=matplotlib.colors.LogNorm(),extend='both')
    cbar2 = fig.colorbar(plot2,ax=ax2,ticks=dum_ticks)
    cbar2.ax.set_ylabel('Rain Rate [mm hr$^{-1}$]')

    ax1.set_title('Native')
    ax2.set_title('Coarse-Grained')
    plt.show()
    plt.close()
        
        
savdir = '/glade/scratch/mckenna/AMIE/amie_post/inst_rain_rates/'
f = open(savdir+'baseline_1km_CG_inst_rain_rates.p','wb')
pickle.dump(rain_rate_deg_dict,f)
f.close()

print('Done')


    
