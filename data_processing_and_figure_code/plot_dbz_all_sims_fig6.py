#==========================================================================
# Title: plot_dbz_all_sims_fig6.py
# Author: McKenna W. Stanford
# Utility: Plots reflectivity of all simulations at 2.5 km AGL.
# Produces Fig. 6 of manuscript.
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
#.use('Agg')
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
import matplotlib.ticker as mticker

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
wrfx0 = np.arange(-148.5,151.5,3)
wrfy0 = np.arange(-148.5,151.5,3)
tmpmeshx,tmpmeshy = np.meshgrid(wrfx0,wrfy0)
#calculate radial distance from origin
radial_dist = np.sqrt(tmpmeshx**2. + tmpmeshy**2.)


#======================================
# read in spol reflectivity
# --not the coarse grained one for now
#======================================
path = '/glade/scratch/mckenna/AMIE/amie_post/spol/'

f = open(path+'spol_coarse_grained.p','rb')
spol_deg = pickle.load(f)['1km']
f.close()

sdbz = spol_deg['dbz']
stime = spol_deg['time']
slon = spol_deg['lon']
slat = spol_deg['lat']
#only grab 2.5 km reflectivity
sdbz = sdbz[:,:,:,4]
stime = stime[60]
sdbz = sdbz[60,:,:]

#======================================
# read in 1 km run
#======================================
path = '/glade/scratch/mckenna/AMIE/baseline_1km/'
files = sorted(glob.glob(path+'wrfout*d02*'))
tmpfile = files[84]

ncfile = Dataset(tmpfile)
lat_1km = wrf.getvar(ncfile,'lat',meta=False).data.T
lon_1km = wrf.getvar(ncfile,'lon',meta=False).data.T
dbz_1km = wrf.getvar(ncfile,'REFL_10CM',meta=False)
z_1km = wrf.getvar(ncfile,'height',meta=False)
ncfile.close()

tmpdbz = wrf.interplevel(dbz_1km,z_1km,2500.)
dbz_1km_fin = tmpdbz.data.T

        
#======================================
# Read in lats/lons from 3-km runs
# and defin lists with parameters
# for subplots
#======================================

sims = ['baseline',\
       'no_mixing',\
       '4x',\
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
       ]
titles = ['BASELINE',\
       'NO\_MIXING',\
       '4X',\
       'STOCH1\_SHORT',\
       'STOCH2\_SHORT',\
       'STOCH3\_SHORT',\
       'STOCH4\_SHORT',\
       'STOCH5\_SHORT',\
       'STOCH1\_LONG',\
       'STOCH2\_LONG',\
       'STOCH3\_LONG',\
       'STOCH4\_LONG',\
       'STOCH5\_LONG',\
       ]

labs = ['(c)','(d)','(e)','(f)','(g)','(h)',\
       '(i)','(j)','(k)','(l)','(m)','(n)','(o)']

path = '/glade/scratch/mckenna/AMIE/baseline/'
files = sorted(glob.glob(path+'wrfout*d01*'))
ncfile = Dataset(files[0])
lon_3km = wrf.getvar(ncfile,'lon',meta=False).data.T
lat_3km = wrf.getvar(ncfile,'lat',meta=False).data.T
ncfile.close()

#======================================
# Start Figure
#======================================

fig = plt.figure(figsize=(12,7))    
axnums = [3,4,5,6,7,8,9,10,11,12,13,14,15]
levs = np.arange(0,62.5,2.5)    
Fontsize=14
fac = 0.75
fac2 = 1.
fac3 = 1.9
dumx = -0.275
dumy = 1.125

dum_lats = np.array([-2,0,2])
dum_lons = np.array([71,73,75])
#-------------------
# SPOL
#-------------------
ax = fig.add_subplot(3,5,1,projection=ccrs.PlateCarree()) 
ax.add_feature(cfeature.OCEAN)
ax.coastlines('10m')
ax.set_xlim(np.min(lon_3km[200-100:200+100,200-100:200+100]),np.max(lon_3km[200-100:200+100,200-100:200+100]))
ax.set_ylim(np.min(lat_3km[200-100:200+100,200-100:200+100]),np.max(lat_3km[200-100:200+100,200-100:200+100]))

g1 = ax.gridlines(color='grey',lw=0.2,\
                  linestyle='solid')
g1.xlabels_bottom = True
g1.ylabels_left = True
g1.xlabel_style = {'size':Fontsize*fac}
g1.ylabel_style = {'size':Fontsize*fac}
g1.xlocator = mticker.FixedLocator(dum_lons)
g1.ylocator = mticker.FixedLocator(dum_lats)
g1.xformatter = LONGITUDE_FORMATTER
g1.yformatter = LATITUDE_FORMATTER

ax.set_aspect('auto')
ax.set_title('SPOL',fontsize=Fontsize*fac2,ha='center')
ax.text(dumx,dumy,'(a)',transform=ax.transAxes,\
           fontsize=Fontsize*fac3,fontweight='bold')
ax.xaxis.set_ticks(dum_lons)
ax.yaxis.set_ticks(dum_lats)

start_id = 200
diff = 50
# Plot range ring
ax.contour(lon_3km[start_id-diff:start_id+diff,\
                start_id-diff:start_id+diff],\
           lat_3km[start_id-diff:start_id+diff,\
                start_id-diff:start_id+diff],radial_dist,\
          levels=[148,151],colors='black',linewidths=2)

ax.contourf(slon,slat,sdbz,\
            cmap=get_cmap('nipy_spectral'),\
            transform=crs.PlateCarree(),\
            levels=levs,\
            norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256))  




#-------------------
# baseline_1km
#-------------------

ax = fig.add_subplot(3,5,2,projection=ccrs.PlateCarree()) 
ax.add_feature(cfeature.OCEAN)
ax.coastlines('10m')
ax.set_xlim(np.min(lon_3km[200-100:200+100,200-100:200+100]),np.max(lon_3km[200-100:200+100,200-100:200+100]))
ax.set_ylim(np.min(lat_3km[200-100:200+100,200-100:200+100]),np.max(lat_3km[200-100:200+100,200-100:200+100]))

g1 = ax.gridlines(color='grey',lw=0.2,\
                  linestyle='solid')
g1.xlabels_bottom = True
g1.ylabels_left = True
g1.xlabel_style = {'size':Fontsize*fac}
g1.ylabel_style = {'size':Fontsize*fac}
g1.xlocator = mticker.FixedLocator(dum_lons)
g1.ylocator = mticker.FixedLocator(dum_lats)
g1.xformatter = LONGITUDE_FORMATTER
g1.yformatter = LATITUDE_FORMATTER

ax.set_aspect('auto')
ax.set_title('BASELINE\_1KM',fontsize=Fontsize*fac2,ha='center')
ax.text(dumx,dumy,'(b)',transform=ax.transAxes,\
           fontsize=Fontsize*fac3,fontweight='bold')


ax.contour(lon_3km[200-50:200+50,200-50:200+50],lat_3km[200-50:200+50,200-50:200+50],radial_dist,\
          levels=[148,151],colors='black',linewidths=2)

ax.contourf(lon_1km,lat_1km,dbz_1km_fin,\
            cmap=get_cmap('nipy_spectral'),\
            transform=crs.PlateCarree(),\
            levels=levs,\
            norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256))  


dum = 0
for ax in range(len(axnums)):
    sim = sims[dum]
    
    print('Simulation: {}'.format(sims[dum]))

    ax = fig.add_subplot(3,5,axnums[dum],projection=ccrs.PlateCarree()) 
    ax.add_feature(cfeature.OCEAN)
    ax.coastlines('10m')
    ax.set_xlim(np.min(lon_3km[200-100:200+100,200-100:200+100]),np.max(lon_3km[200-100:200+100,200-100:200+100]))
    ax.set_ylim(np.min(lat_3km[200-100:200+100,200-100:200+100]),np.max(lat_3km[200-100:200+100,200-100:200+100]))

    g1 = ax.gridlines(color='grey',lw=0.2,linestyle='solid')
    g1.xlabels_bottom = True
    g1.ylabels_left = True
    g1.xlabel_style = {'size':Fontsize*fac}
    g1.ylabel_style = {'size':Fontsize*fac}
    ax.set_aspect('auto')
    g1.xlocator = mticker.FixedLocator(dum_lons)
    g1.ylocator = mticker.FixedLocator(dum_lats)
    g1.xformatter = LONGITUDE_FORMATTER
    g1.yformatter = LATITUDE_FORMATTER
    
    
    ax.set_title(titles[dum],fontsize=Fontsize*fac2,ha='center')
    ax.text(dumx,dumy,labs[dum],transform=ax.transAxes,\
           fontsize=Fontsize*fac3,fontweight='bold')    
    
    
    # plot range ring
    ax.contour(lon_3km[200-50:200+50,200-50:200+50],lat_3km[200-50:200+50,200-50:200+50],radial_dist,\
              levels=[148,151],colors='black',linewidths=2)

    # Read in file
    tmpfile = sorted(glob.glob('/glade/scratch/mckenna/AMIE/'+sim+'/wrfout_d01*'))[108]


    ncfile = Dataset(tmpfile)
    dbz_3km = wrf.getvar(ncfile,'REFL_10CM',meta=False)
    z_3km = wrf.getvar(ncfile,'height',meta=False)
    ncfile.close()

    tmpdbz = wrf.interplevel(dbz_3km,z_3km,2500.)
    dbz_3km_fin = tmpdbz.data.T


    ax.contourf(lon_3km,lat_3km,dbz_3km_fin,\
            cmap=get_cmap('nipy_spectral'),\
            transform=crs.PlateCarree(),\
            levels=levs,\
            norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256))               

    dum+=1
    

cbar_ax = fig.add_axes([0.13,0.05,0.75,0.03])
# Set color mappable
range_min = 0
range_max = 60
cmap = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256),\
cmap = plt.get_cmap('nipy_spectral'))    

cmap.set_array([])   
cbar = fig.colorbar(cmap,orientation='horizontal',cax=cbar_ax)
cbar.ax.tick_params(labelsize=Fontsize*1.4)
cbar.set_label('Rayleigh Reflectivity [dBZ]',fontsize=Fontsize*1.4)

plt.subplots_adjust(hspace=0.4,wspace = 0.25, top=0.915)

plt.suptitle('0300 UTC, 8 Dec. 2.5 km AGL Reflectivity',fontsize=Fontsize*1.8,\
            ha='center',va='bottom')


save_path = '/glade/u/home/mckenna/figures/amie_paper/'
outfile = 'fig_06.png'
plt.savefig(save_path+outfile,dpi=200,bbox_inches='tight')
#plt.show()
plt.close()
print('done')
    
    









