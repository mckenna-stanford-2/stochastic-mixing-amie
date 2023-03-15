#==========================================================================
# Title: plot_amie_wrf_dbz.py
# Author: McKenna W. Stanford
# Utility: Plots reflectivity cross sections of observed and simulated
# radar reflectivity.
# Produces Fig. 4 of the manuscripts.
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
import wrf
from netCDF4 import Dataset
import cartopy.crs as crs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.offsetbox import AnchoredText
import datetime
import wrf
from netCDF4 import Dataset
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
#------------------------------------------
# Parameters
#------------------------------------------
plt.rc('text', usetex=True)

#------------------------------------------
# Functions
#------------------------------------------
# make function that finds the nearest
# element in array to the given value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx


# Define SPOL range
wrfx0 = np.arange(-149.5,150.5,1)
wrfy0 = np.arange(-149.5,150.5,1)
tmpmeshx,tmpmeshy = np.meshgrid(wrfx0,wrfy0)
#calculate radial distance from origin
radial_dist = np.sqrt(tmpmeshx**2. + tmpmeshy**2.)


iread = 0
if iread == 1.:

    #=================================
    # SPOL
    #=================================
    path = '/glade/scratch/mckenna/AMIE/spol/refl/'
    spol_files = np.array(sorted(glob.glob(path+'*.nc')))

    #sids = [68,92,116,140]
    sids = [80,92,116,140]
    spol_files = spol_files[sids]

    sdbz = []
    for file in spol_files:
        print(file)
        ncfile = xarray.open_dataset(file)
        slon = np.array(ncfile['lon0'].copy())
        slat = np.array(ncfile['lat0'].copy())
        tmp_sdbz = np.array(ncfile['REFL'].copy())
        sz = np.array(ncfile['z0'].copy())
        ncfile.close()
        slat = np.transpose(slat)
        slon = np.transpose(slon)
        tmp_sdbz = np.squeeze(tmp_sdbz)
        tmp_sdbz = np.transpose(tmp_sdbz)

        tmpz,tmpzid = find_nearest(sz,2.5)
        tmp_sdbz = tmp_sdbz[:,:,tmpzid]
        sdbz.append(tmp_sdbz)
    sdbz = np.array(sdbz)


    #=================================
    # 3-km BASELINE run
    #=================================
    path1 = '/glade/scratch/mckenna/AMIE/baseline/'
    files1 = glob.glob(path1+'wrfout_d01*')
    files1 = sorted(files1)
    #ids1 = [72,96,120,144]
    ids1 = [84,96,120,144]
    files1 = np.array(files1)
    files1 = files1[ids1]
    #print(aaaa)

    dbz_3km = []
    for file in files1:
        print(file)

        ncfile = Dataset(file)
        z = wrf.getvar(ncfile,'height',meta=False)
        lat3 = wrf.getvar(ncfile,'lat',meta=False)
        lon3 = wrf.getvar(ncfile,'lon',meta=False)
        tmpdbz = wrf.getvar(ncfile,'REFL_10CM',meta=False)
        ncfile.close()

        tmpze = tmpdbz.copy()
        dumid = np.where(tmpdbz == -99.)
        tmpze[dumid] = 0.
        dumid = np.where(tmpdbz > -99.)
        tmpze[dumid] = 10.**(tmpdbz[dumid]/10.)
        tmpze = wrf.interplevel(tmpze,z,2500.)
        tmpdbz = 10.*np.log10(tmpze)
        
        lat3 = lat3.data.T
        lon3 = lon3.data.T
        tmpdbz = tmpdbz.data.T

        dbz_3km.append(tmpdbz)
    dbz_3km = np.array(dbz_3km)
    print('done')
    #plt.contourf(lon3,lat3,dbz_3km[3,:,:],levels=np.arange(0,65,5),cmap='nipy_spectral')

    #=================================
    # 1-km BASELINE run
    #=================================
    path2 = '/glade/scratch/mckenna/AMIE/baseline_1km/'
    files2 = glob.glob(path2+'wrfout_d02*:00')
    files2 = sorted(files2)
    ids2 = [60,72,96,120]
    #ids2 = [60,72,96,115]
    files2 = np.array(files2)
    files2 = files2[ids2]

    dbz_1km = []
    for file in files2:
        print(file)


        ncfile = Dataset(file)
        z = wrf.getvar(ncfile,'height',meta=False)
        lat1 = wrf.getvar(ncfile,'lat',meta=False)
        lon1 = wrf.getvar(ncfile,'lon',meta=False)
        tmpdbz = wrf.getvar(ncfile,'REFL_10CM',meta=False)
        ncfile.close()
        
        tmpze = tmpdbz.copy()
        dumid = np.where(tmpdbz == -99.)
        tmpze[dumid] = 0.
        dumid = np.where(tmpdbz > -99.)
        tmpze[dumid] = 10.**(tmpdbz[dumid]/10.)
        tmpze = wrf.interplevel(tmpze,z,2500.)
        tmpdbz = 10.*np.log10(tmpze)
        
        lat1 = lat1.data.T
        lon1 = lon1.data.T
        tmpdbz = tmpdbz.data.T
        
        dbz_1km.append(tmpdbz)
    dbz_1km = np.array(dbz_1km)
    #plt.contourf(lon1,lat1,dbz_1km[3,:,:],levels=np.arange(0,65,5),cmap='nipy_spectral')










#====================================================================
#====================================================================
#====================================================================
# Plots
#====================================================================
#====================================================================
#====================================================================


Fontsize=16

#============================
# Start Figure
#============================
fig = plt.figure(figsize=(12,16))
ax1 = fig.add_subplot(431,projection=ccrs.PlateCarree())
ax2 = fig.add_subplot(432,projection=ccrs.PlateCarree())
ax3 = fig.add_subplot(433,projection=ccrs.PlateCarree())
ax4 = fig.add_subplot(434,projection=ccrs.PlateCarree())
ax5 = fig.add_subplot(435,projection=ccrs.PlateCarree())
ax6 = fig.add_subplot(436,projection=ccrs.PlateCarree())
ax7 = fig.add_subplot(437,projection=ccrs.PlateCarree())
ax8 = fig.add_subplot(438,projection=ccrs.PlateCarree())
ax9 = fig.add_subplot(439,projection=ccrs.PlateCarree())
ax10 = fig.add_subplot(4,3,10,projection=ccrs.PlateCarree())
ax11 = fig.add_subplot(4,3,11,projection=ccrs.PlateCarree())
ax12 = fig.add_subplot(4,3,12,projection=ccrs.PlateCarree())

axlist = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12]
levs = np.arange(0,62.5,2.5)

titles = ['7 Dec., 21 UTC',\
         '7 Dec., 21 UTC',\
         '7 Dec., 21 UTC',\
         '8 Dec., 00 UTC',\
         '8 Dec., 00 UTC',\
         '8 Dec., 00 UTC',\
         '8 Dec., 06 UTC',\
         '8 Dec., 06 UTC',\
         '8 Dec., 06 UTC',\
         '8 Dec., 12 UTC',\
         '8 Dec., 12 UTC',\
         '8 Dec., 12 UTC',\
        ]
dum = 0
for ax in axlist:
    ax.add_feature(cfeature.OCEAN)
    ax.coastlines('10m')
    ax.set_xlim(np.min(slon)-0.3,np.max(slon)+0.3)
    ax.set_ylim(np.min(slat)-0.3,np.max(slat)+0.3)
    g1 = ax.gridlines(color='black',linestyle='solid',lw=1.)
    g1.xlabels_bottom = True
    g1.ylabels_left = True
    g1.xlabel_style = {'size':Fontsize}
    g1.ylabel_style = {'size':Fontsize}
    ax.set_aspect('auto')
     
    dum = dum + 1
    

# SPOL plots
axlist0 = [ax1,ax4,ax7,ax10]
dum = 0
for ax in axlist0:
    ax.contourf(slon,slat,sdbz[dum,:,:],\
        cmap=get_cmap('nipy_spectral'),\
        transform=crs.PlateCarree(),\
        levels=levs,\
        norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256))
    ax.contour(slon,slat,\
               radial_dist,\
               levels=[149,151],colors='black',linewidths=5)    
    dum+=1
    
# 1-km plots
axlist1 = [ax2,ax5,ax8,ax11]
dum = 0
for ax in axlist1:
    ax.contourf(lon1,lat1,dbz_1km[dum,:,:],\
        cmap=get_cmap('nipy_spectral'),\
        transform=crs.PlateCarree(),\
        levels=levs,\
        norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256))
    ax.contour(slon,slat,\
               radial_dist,\
               levels=[149,151],colors='black',linewidths=5)
    dum+=1
    
# 3-km plots
axlist2 = [ax3,ax6,ax9,ax12]
dum = 0
for ax in axlist2:
    ax.contourf(lon3,lat3,dbz_3km[dum,:,:],\
        cmap=get_cmap('nipy_spectral'),\
        transform=crs.PlateCarree(),\
        levels=levs,\
        norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256))
    ax.contour(slon,slat,\
               radial_dist,\
               levels=[149,151],colors='black',linewidths=5)    
    dum+=1
        
    
    
cbar_ax = fig.add_axes([0.13,0.07,0.75,0.03])   
# Set color mappable
range_min = 0
range_max = 60
cmap = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256),\
    cmap = plt.get_cmap('nipy_spectral'))    
cmap.set_array([])   
    
cbar = fig.colorbar(cmap,orientation='horizontal',cax=cbar_ax)
cbar.ax.tick_params(labelsize=Fontsize*1.5)
cbar.set_label('2.5-km Rayleigh Reflectivity [dBZ]',fontsize=Fontsize*1.5)

plt.text(0.5,1.1,'SPOL',fontsize=Fontsize*2,fontweight='bold',\
         #bbox=dict(facecolor='grey', alpha=0.25),transform=ax1.transAxes,\
         transform=ax1.transAxes,\
        ha='center')
plt.text(0.5,1.1,'$\\Delta_{h}$ = 1 km',fontsize=Fontsize*2,fontweight='bold',\
         #bbox=dict(facecolor='grey', alpha=0.25),transform=ax2.transAxes,\
         transform=ax2.transAxes,\
        ha='center')
plt.text(0.5,1.1,'$\\Delta_{h}$ = 3 km',fontsize=Fontsize*2,fontweight='bold',\
         #bbox=dict(facecolor='grey', alpha=0.25),transform=ax3.transAxes,\
         transform=ax3.transAxes,\
        ha='center')

plt.text(-0.35,0.5,'7 Dec., 21Z',fontsize=Fontsize*2,fontweight='bold',\
         #bbox=dict(facecolor='deepskyblue', alpha=0.25),transform=ax1.transAxes,\
         transform=ax1.transAxes,\
         rotation=90,va='center')
plt.text(-0.35,0.5,'8 Dec., 00Z',fontsize=Fontsize*2,fontweight='bold',\
         #bbox=dict(facecolor='deepskyblue', alpha=0.25),transform=ax4.transAxes,\
         transform=ax4.transAxes,\
         rotation=90,va='center')
plt.text(-0.35,0.5,'8 Dec., 06Z',fontsize=Fontsize*2,fontweight='bold',\
         #bbox=dict(facecolor='deepskyblue', alpha=0.25),transform=ax7.transAxes,\
         transform=ax7.transAxes,\
         rotation=90,va='center')
plt.text(-0.35,0.5,'8 Dec., 12Z',fontsize=Fontsize*2,fontweight='bold',\
         #bbox=dict(facecolor='deepskyblue', alpha=0.25),transform=ax10.transAxes,\
         transform=ax10.transAxes,\
         rotation=90,va='center')



labs = ['(a)','(b)','(c)','(d)','(e)',\
       '(f)','(g)','(h)','(i)','(j)','(k)','(l)']

dum = 0
for ax in axlist:
    ax.text(0.025,0.85,labs[dum],transform=ax.transAxes,fontweight='bold',\
            bbox=dict(facecolor='lightgrey', alpha=0.75),fontsize=Fontsize*2.)
    dum+=1

plt.subplots_adjust(wspace=0.25,hspace=0.25)    
    

outfile = 'fig_04.png'
save_path = '/glade/u/home/mckenna/figures/amie_paper/'
plt.savefig(save_path+outfile,dpi=260,bbox_inches='tight')
#plt.show()
plt.close()

#===============================================================
#===============================================================
#===============================================================
#===============================================================
#===============================================================



