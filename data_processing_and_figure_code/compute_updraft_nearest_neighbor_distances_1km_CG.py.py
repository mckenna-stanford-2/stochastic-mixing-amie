#==========================================================================
# Title: compute_updraft_nearest_neighbor_distances_1km_CG.py
# Author: McKenna W. Stanford
# Utility: Finds updraft objects and computes the distance to the
# nearest in space updraft object. 
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
import pandas as pd
from scipy.ndimage import label
import cartopy.crs as ccrs
from matplotlib.cm import get_cmap
from matplotlib.lines import Line2D
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


path = '/glade/scratch/mckenna/AMIE/'

files = sorted(glob.glob(path+'baseline_1km/3km/'+'/wrfout*'))
nt = len(files)

nnd_all = []
#========================================================
# Loop through times
#========================================================
#for tt in range(0,5):
for tt in range(nt):
    print('Time Step: {}/{}'.format(tt+1,nt))
    ncfile = xarray.open_dataset(files[tt])
    w = ncfile['w'].values
    z = ncfile['z'].values
    twc = ncfile['twc'].values
    lon = ncfile['lon'].values
    lat = ncfile['lat'].values

    nx = len(w[:,0,0])
    ny = len(w[0,:,0])
    nz = len(w[0,0,:])

    avg_z = np.nanmean(z,axis=(0,1))        
    z_25km,id_25km = find_nearest(avg_z,2500.)
    w = w[:,:,id_25km]
    twc = twc[:,:,id_25km]
    up_mask = np.zeros((nx,ny))

    where_gt_1 = np.where((w > 1.) & (twc > 0.1))
    up_mask[where_gt_1] = 1

    regions = label(up_mask)
    label_num = regions[1]
    regions = regions[0]


    #--------------------------------------------
    # Loop through updraft objects and calculate
    # the center of mass
    #--------------------------------------------
    labeled_updrafts = np.zeros((nx,ny))
    center_of_mass = np.zeros((label_num,2))
    center_of_mass_id = np.zeros((label_num,2))

    dum = 0
    for tmplabel in range(label_num):
        tmpid = np.where(regions == (tmplabel + 1))
        labeled_updrafts[tmpid] = tmplabel + 1
        object_lon = lon[tmpid]
        object_lat = lat[tmpid]  
        object_w = w[tmpid]

        xcm = np.sum(object_w*object_lon)/np.sum(object_w)
        ycm = np.sum(object_w*object_lat)/np.sum(object_w)
        nearest_xcm, nearest_xid = find_nearest(object_lon,xcm)
        nearest_ycm, nearest_yid = find_nearest(object_lat,ycm)   
        center_of_mass[dum,:] = np.array([nearest_xcm,nearest_ycm])

        tmp_lon = object_lon[nearest_xid]
        tmp_lat = object_lat[nearest_yid]

        newid = np.where((lon == tmp_lon) & (lat == tmp_lat))
        center_of_mass_id[dum,:] = np.array([newid[0][0],newid[1][0]])


        #------------------
        # plot
        #------------------
        if False:
            fig = plt.figure(figsize=(8,6))
            ax1 = fig.add_subplot(111)
            tmpplot=ax1.contourf(lon,lat,w,levels=np.arange(-5,5.1,0.1),\
                cmap='seismic',extend='both')
            ax1.scatter(lon[tmpid],lat[tmpid],200,color='green',marker='o')

            xloc_mean = int(np.round(np.mean(tmpid[0])))
            yloc_mean = int(np.round(np.mean(tmpid[1])))
            xloc_min = int(np.round(np.min(tmpid[0])))
            xloc_max = int(np.round(np.max(tmpid[0])))
            yloc_min = int(np.round(np.min(tmpid[1])))
            yloc_max = int(np.round(np.max(tmpid[1])))  

            ax1.set_xlim(lon[xloc_min,yloc_min]-0.05,lon[xloc_max,yloc_max]+0.05)
            ax1.set_ylim(lat[xloc_min,yloc_min]-0.05,lat[xloc_max,yloc_max]+0.05)
            cbar = plt.colorbar(tmpplot,ticks=np.arange(-5,6,1))
            cbar.ax.set_ylabel('$w$ [m s$^{-1}$]',fontsize=20)
            cbar.ax.tick_params(labelsize=15)
            plt.tick_params(labelsize=15)

             # plot "x" where center of mass is
            ax1.plot(center_of_mass[dum][0],center_of_mass[dum][1],'x',color='blue',\
                     markersize=20,markeredgewidth=7.5,label='Center of Mass')           

            plt.show()

        dum+=1

    #--------------------------------------------
    # Loop through updraft objects and calculate
    # the nearest neighbor distance
    #--------------------------------------------
    nnd = []
    # look through updrafts
    for up in range(label_num):
        # Get center of mass and the id of it
        tmp_center_of_mass = center_of_mass[up,:]
        tmp_center_of_mass_id = center_of_mass_id[up,:]
        # calcualte distance between the id of the current
        # center of mass and 
        xdiff = abs(tmp_center_of_mass_id[0]-center_of_mass_id[:,0])
        ydiff = abs(tmp_center_of_mass_id[1]-center_of_mass_id[:,1])
        # calculate radial distance
        rad_diff = np.sqrt(xdiff**2. + ydiff**2.)
        # grab nearest neighbor distance that is
        # greater than 0 (because zero is the current)
        # point
        nn = np.min(rad_diff[rad_diff > 0.])
        nnd.append(nn)

    # Array of  distance between each updraft object
    # and the closest object to it.
    nnd = np.array(nnd)
    # Append to array
    nnd_all.append(nnd)        

    #--------------------------------------------        
    # plot empirical CDF of nearest neighbor
    # distances
    #--------------------------------------------        
    iplot = False
    if iplot:
        fig = plt.figure(figsize=(16,8))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        Fontsize=22
        ax1.grid(which='both',ls='dotted',lw=1,c='grey')
        ax2.grid(which='both',ls='dotted',lw=1,c='grey')
        x = np.sort(nnd)
        y = np.arange(1,len(x)+1)/len(x)
        # Plot CDF of nearest neighbor unit distances
        ax1.plot(x,y,label='NNCDF')

        r = np.sort(x)
        L = 400
        N = len(x)
        lam = N/(L**2.)
        # Plot Random Distribution
        yran = 1 - np.exp(-1.*lam*np.pi*(r**2.))
        ax1.plot(r,yran,label='NNCDF$_{ran}$')

        # Plot nearest neighbor CDF as a function
        # of completely random CDF
        ax2.plot(yran,y)
        # Plot relationship for completely random distribution
        ax2.plot(np.arange(0,1.1,0.1),np.arange(0,1.1,0.1),c='k')

        ax1.legend(fontsize=Fontsize)
        ax1.set_ylabel('NNCDF',fontsize=Fontsize)
        ax1.set_xlabel('Nearest Neighbor UnitDistance',fontsize=Fontsize)
        ax2.set_ylabel('NNCDF',fontsize=Fontsize)
        ax2.set_xlabel('NNCDF$_{ran}$',fontsize=Fontsize)
        ax1.tick_params(labelsize=Fontsize)
        ax2.tick_params(labelsize=Fontsize)

        plt.show()
        plt.close()

# Now we have the nearest neighbor distances as a list
# for each time period. Let's create a dictionary
# indexed by a string representing the time stamp
nnd_dict = {}
for ii in range(len(nnd_all)):
    nnd_dict[str(ii+1)] = nnd_all[ii]

for key,val in nnd_dict.items():
    print(key,np.shape(val))

#--------------------------------------------        
# Calculate nearest neighbor distance CDFs
# and intgrate them
#--------------------------------------------  

# Save to dictionary for each simulation
print('Writing to dictionary...')
savdir = '/glade/scratch/mckenna/AMIE/amie_post/updraft_nearest_neighbor/'
f = open(savdir+'baseline_1km_CG_updraft_nearest_neighbor.p','wb')
pickle.dump(nnd_dict,f)
f.close()
print('Completed writing to dictionary. Done.')
    
    
    
print(aaaaa)
#------------------------------------
# Plots
#------------------------------------
if False:
    def cdf(data):

        data_size=len(data)

        # Set bins edges
        data_set=sorted(set(data))
        #bins=np.append(data_set, data_set[-1]+1)
        bins = np.arange(0,402,2)

        # Use the histogram function to bin the data
        counts, bin_edges = np.histogram(data, bins=bins, density=False)

        counts=counts.astype(float)/data_size

        # Find the cdf
        cdf = np.cumsum(counts)

        return cdf, bin_edges


    nt = len(time)

    for tt in range(nt):


        plt.rc('text',usetex=True)
        fig = plt.figure(figsize=(10,4))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.grid(which='both',ls='dotted',lw=1,c='grey')
        ax2.grid(which='both',ls='dotted',lw=1,c='grey')
        ax1.set_ylabel('NNCDF',fontsize=15)
        ax1.set_xlabel('Nearest Neighbor Distance [km]',fontsize=15)
        ax2.set_ylabel('NNCDF',fontsize=15)
        ax2.set_xlabel('NNCDF$_{ran}$',fontsize=15) 
        ax2.plot(np.arange(0,1.1,0.1),np.arange(0,1.1,0.1),c='grey',lw=3)
        ax1.tick_params(labelsize=15)
        ax2.tick_params(labelsize=15)
        ax1.set_xlim(0,150)


        dum = 0
        colors=['black','red','blue']
        for key,val in nnd_dict.items():
            print(key,np.shape(val[tt]))

            x = val[tt]*3.
            #x = np.sort(val[tt])*3.
            y,bin_edges = cdf(x)
            #y = np.arange(1,len(x)+1)/len(x)
            #ax1.plot(x,y,color=colors[dum],lw=2)
            mid_bins = []
            for ii in range(len(bin_edges)-1):
                mid_bins.append(0.5*(bin_edges[ii+1]+bin_edges[ii]))

            #ax1.plot(bin_edges[0:-1],y,color=colors[dum],lw=2)
            ax1.plot(mid_bins,y,color=colors[dum],lw=2)

            #r = np.sort(bin_edges[0:-1])
            r = np.sort(mid_bins)
            L = 400*3.
            N = len(x)
            lam = N/(L**2.)
            yran = 1 - np.exp(-1.*lam*np.pi*(r**2.))

            ax2.plot(yran,y,color=colors[dum],lw=2)

            dum+=1

        ax1.plot(r,yran,color='grey',lw=2)


        custom_lines_1 = [Line2D([0],[0],color = 'black',lw=3),\
                        Line2D([0],[0],color='blue',ls='solid',lw=3),\
                        Line2D([0],[0],color='red',ls='solid',lw=3),\
                        Line2D([0],[0],color='grey',ls='solid',lw=3),\
                       ]

        ax1.legend(custom_lines_1,['BASELINE','4X','NO\_MIXING','Theoretical'],fontsize=12,\
                        loc='lower right',ncol=1)  


        custom_lines_2 = [Line2D([0],[0],color = 'black',lw=3),\
                        Line2D([0],[0],color='blue',ls='solid',lw=3),\
                        Line2D([0],[0],color='red',ls='solid',lw=3),\
                       ]
        ax2.legend(custom_lines_2,['BASELINE','4X','NO\_MIXING'],fontsize=12,\
                        loc='lower right',ncol=1)   
        plt.subplots_adjust(wspace=0.3)
        plt.suptitle('time = {} UTC'.format(str(time[tt])),fontsize=25)

        savpath = '/glade/u/home/mckenna/figures/nncdfs/'
        timestr = np.around(time2[tt],4)
        outfile = 'nncdfs_{}.png'.format(timestr)
        plt.savefig(savpath+outfile,dpi=160,bbox_inches='tight')
        #plt.show()    

        #print(aaaaa)

    
import pickle
savdir = '/glade/scratch/mckenna/AMIE/amie_post/'
#f = open(savdir+'nnd_dict','wb')
f = open(savdir+'nnd_dict','wb')
pickle.dump(nnd_dict,f)
f.close()

# test read in
#pkl_file = open(savdir+'nnd_dict','rb')
#new_nnd_dict = pickle.load(pkl_file)
#pkl_file.close()
