#==========================================================================
# Title: compute_spol_conv_strat_regions.py
# Author: McKenna W. Stanford
# Utility: Computes convective and stratiform regions on SPOL data
# following algorithm of Steiner et al. (1995) and saves masks, dbz, and 
# rain rates @ 2.5 km AGL to a dictionary named spol_conv_strat_id.p'
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
import cartopy.crs as ccrs
from matplotlib.cm import get_cmap

# Read in SPOL variables
path = '/glade/scratch/mckenna/AMIE/amie_post/spol/'

f = open(path+'spol_coarse_grained.p','rb')
spol_dict = pickle.load(f)
f.close()

# Print diagnostics
if False:
#if True:
    for key,val in spol_dict.items():
        print(key)
        
        for key2,val2 in spol_dict[key].items():
            print('   ',key2,np.shape(val2))
    

time = spol_dict['3km']['time']
lon = spol_dict['3km']['lon']
lat = spol_dict['3km']['lat']
nt = len(time)

#------------------------------------------
# Function that finds the nearest
# element in array to the given value
#------------------------------------------

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

#========================================================
#========================================================
object_dict = {}
#========================================================
#========================================================

con1_mask_arr = []
con2_mask_arr = []
con3_mask_arr = []
dbz_2pt5km_arr = []

#========================================================
# Loop through times
#========================================================
for tt in range(nt):
    print('Time step {}/{}:'.format(tt+1,nt))
    rain_rate = spol_dict['3km']['rain_rate'][tt,:,:]
    dbz_2pt5km = spol_dict['3km']['dbz'][tt,:,:,4] # get it at 2.5 km
    z = spol_dict['3km']['z']
    
    dbz_2pt5km_arr.append(dbz_2pt5km)
    
    x0 = np.arange(-148.5,151.5,3)
    y0 = np.arange(-148.5,151.5,3)
    tmpmeshx,tmpmeshy = np.meshgrid(x0,y0)
    #calculate radial distance from origin
    radial_dist = np.sqrt(tmpmeshx**2. + tmpmeshy**2.)

    nx = len(dbz_2pt5km[:,0])
    ny = len(dbz_2pt5km[0,:])
    
    dx_km = 3.
    dy_km = 3.

    tmpid = np.where(~np.isnan(dbz_2pt5km))
    ze_2pt5km = np.zeros((nx,ny))
    ze_2pt5km[tmpid] = 10.**(dbz_2pt5km[tmpid]/10.)  
        
    dbz_con1 = np.zeros((nx,ny))
    dbz_con2 = np.zeros((nx,ny))
    dbz_con3 = np.zeros((nx,ny))
    dbz_background = np.zeros((nx,ny))-999.
    
    
    #------------------------------------------------------
    #------------------------------------------------------
    #------------------------------------------------------
    # STEP 1: Take all cells where dbz > 40
    #------------------------------------------------------
    #------------------------------------------------------
    #------------------------------------------------------
    print('Starting step 1...')
    where_high = np.where(dbz_2pt5km > 40.)
    dbz_con1[where_high] = 1

    regions_40 = label(dbz_con1)
    label_num = regions_40[1]
    regions_40 = regions_40[0]
    print('Completed step 1. # of objects w/ reflectivity > 40 dbz:',label_num)   
            
        
    #------------------------------------------------------
    #------------------------------------------------------
    #------------------------------------------------------
    # STEPS 2 & 3: 
    # Step 2: Take cells within 11-km radius where dbz >
    # locally-dependent threshold
    # Step 3: Laterally dxpand regions within locally-determined
    # radius.
    # Combined since we are looping through horizontal domain.
    #------------------------------------------------------
    #------------------------------------------------------
    #------------------------------------------------------
    print('Starting Steps 2 and 3...')

    bg_thresh = 11.

    dbz_con2 = dbz_con1.copy() # start with Step 1 values
    dbz_con3 = dbz_con1.copy() # start with Step 2 values

    for ix in range(nx):
        for iy in range(ny):
            # find locally-dependent trehsold value
            dbz_bg = dbz_con1.copy()*0.
            ix_min = max(ix-int(bg_thresh/dx_km),0)
            ix_max = min(ix+int(bg_thresh/dx_km),nx-1)
            iy_min = max(iy-int(bg_thresh/dx_km),0)
            iy_max = min(iy+int(bg_thresh/dx_km),ny-1)
            for ix_bg in range(ix_min,(ix_max+1)):
                for iy_bg in range(iy_min,(iy_max+1)):
                    radial_distance = np.sqrt(((dx_km*(ix-ix_bg))**2.)+((dx_km*(iy-iy_bg))**2.))
                    if radial_distance <= bg_thresh:
                        dbz_bg[ix_bg,iy_bg] = 1
            dbz_bg[ix,iy] = 0. # exclude center point itself
            dbz_bg_mean = 0.

            #lt_zero = np.where((dbz_bg > 0) & (ze_2pt5km < 1))
            #lt_zero_size = np.size(lt_zero)
            #if lt_zero_size > 0:
            #    ze_2pt5km[lt_zero] = 1
            where_bg_echo = np.where((dbz_bg > 0) & (ze_2pt5km > 1))
            where_bg_echo_size = np.size(where_bg_echo)
            if where_bg_echo_size > 0:
                ze_bg_mean = np.mean(ze_2pt5km[where_bg_echo])
                dbz_bg_mean = 10.*np.log10(ze_bg_mean)
                dbz_background[ix,iy] = dbz_bg_mean
            else:
                dbz_background[ix,iy] = np.nan
                dbz_bg_mean = np.nan

            if ~np.isnan(dbz_bg_mean):
                if (dbz_bg_mean <= 0.): # different from steiner
                    del_dbz_thresh = 10.
                if ((dbz_bg_mean <= 42.43) & (dbz_bg_mean > 0)):
                    del_dbz_thresh = 10. - ((dbz_bg_mean**2.)/180.)
                if (dbz_bg_mean > 42.43):
                    del_dbz_thresh = 0.
                dbz_thresh  = dbz_bg_mean + del_dbz_thresh # include cells in radius that exceed this dbz


                # apply locally-dependent threshold value
                if ((dbz_2pt5km[ix,iy] > dbz_thresh) & (dbz_2pt5km[ix,iy] >= 0.)):
                    dbz_con2[ix,iy] = 1 # add new values to first round of values
                    dbz_con3[ix,iy] = 1 # add new values to first round of values

                if (dbz_con2[ix,iy] == 1):

                    rad_thresh = 0.
                    if dbz_bg_mean < 20.:
                        rad_thresh = 1.
                    if (dbz_bg_mean >= 20.) & (dbz_bg_mean < 25):
                        rad_thresh = 2
                    if (dbz_bg_mean >= 25.) & (dbz_bg_mean < 30):
                        rad_thresh = 3.
                    if (dbz_bg_mean >= 30) & (dbz_bg_mean < 35.):
                        rad_thresh = 4
                    if (dbz_bg_mean >= 35) & (dbz_bg_mean < 40.):
                        rad_thresh = 5                    
                    if (dbz_bg_mean >= 40.):
                        rad_thresh = 6.                          

                    # find values within local radius
                    if (rad_thresh > 0.):
                        dbz_bg = dbz_con2.copy()*0.
                        ix_min = max(ix-int(6/dx_km),0)
                        ix_max = min(ix+int(6/dx_km),nx-1)
                        iy_min = max(iy-int(6/dx_km),0)
                        iy_max = min(iy+int(6/dx_km),ny-1)
                        for ix_bg in range(ix_min,(ix_max+1)):
                            for iy_bg in range(iy_min,(iy_max+1)):
                                radial_distance = np.sqrt(((dx_km*(ix-ix_bg))**2.)+((dx_km*(iy-iy_bg))**2.))
                                if radial_distance <= rad_thresh:
                                    dbz_bg[ix_bg,iy_bg] = 1 
                        dbz_bg[ix,iy] = 0. # exclude center point itself
                        where_high = np.where(dbz_bg > 0)
                        where_high_size = np.size(where_high)
                        if where_high_size > 0:
                            dbz_con3[where_high] = 1
                        where_not = np.where((dbz_bg > 0) & (ze_2pt5km <= 1.))
                        where_not_size = np.size(where_not)
                        if where_not_size > 0.:
                            dbz_con3[where_not] = 0.
            else:
                continue



    regions_bg = label(dbz_con2)
    label_num_bg = regions_bg[1]
    regions_bg = regions_bg[0]
    print('Completed step 2.')
    print('# of objects meeting background reflectivity threshold requirements:',label_num_bg)  


    lat_expand_regions = label(dbz_con3)
    lat_expand_label_num = lat_expand_regions[1]
    lat_expand_regions = lat_expand_regions[0]    
    print('Completed step 3.') 
    print('# of objects laterally expanding convective regions:',lat_expand_label_num)  


    #---------------------------------------------
    # Identify stratiform regions, as well
    #---------------------------------------------
    where_strat = np.where((dbz_2pt5km > 0) & (dbz_con3 != 1))
    size_where_strat = np.size(where_strat)
    if size_where_strat > 0:
        dbz_con3[where_strat] = 2

    con1_mask_arr.append(dbz_con1)
    con2_mask_arr.append(dbz_con2)
    con3_mask_arr.append(dbz_con3)

    iplot = False
    if iplot:
        tmp_conv = dbz_con3*0.
        tmp_strat = dbz_con3*0.
        tmp_strat[dbz_con3 == 2] = 1
        tmp_conv[dbz_con3 == 1] = 1
        fig = plt.figure(figsize=(18,6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        Fontsize=14
        ax1.contourf(lon,lat,dbz_2pt5km,cmap='nipy_spectral',levels=np.arange(0,65,5))
        ax2.contourf(lon,lat,tmp_strat,colors='blue',levels=[0.5,1.5],label='Strat.')
        tmpplot = ax2.contourf(lon,lat,tmp_conv,colors='red',levels=[0.5,1.5],label='Conv.')
        lgnd = ax2.legend(fontsize=Fontsize,loc='upper right')
        #ax2.scatter(lon,lat,tmp_conv,c='red',label='Conv')
        #ax2.scatter(lon,lat,tmp_strat,c='blue',label='strat')
        plt.show()
        plt.close()

con1_mask_arr = np.array(con1_mask_arr)
con2_mask_arr = np.array(con2_mask_arr)
con3_mask_arr = np.array(con3_mask_arr)
dbz_2pt5km_arr = np.array(dbz_2pt5km_arr)

object_dict['con1_mask'] = con1_mask_arr
object_dict['con2_mask'] = con2_mask_arr
object_dict['con3_mask'] = con3_mask_arr
object_dict['rain_rate'] = spol_dict['3km']['rain_rate']
object_dict['radial_dist'] = radial_dist
object_dict['dbz'] = dbz_2pt5km_arr
object_dict['lon'] = lon
object_dict['lat'] = lat
object_dict['time'] = time

# Print diagnostics
if True:
#if False:
    for key,val in object_dict.items():
        print(key,np.shape(val),np.nanmax(val),np.nanmin(val))
          

print('Writing to dictionary...')
savdir = '/glade/scratch/mckenna/AMIE/amie_post/conv_strat_id/'
f = open(savdir+'spol_conv_strat_id.p','wb')
pickle.dump(object_dict,f)
f.close()    
print('Done.')    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
