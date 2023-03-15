#==========================================================================
# Title: plot_amie_domain_fig 1.py
# Utility: Plots the AMIE domain used to produce Fig. 01. 
# Author: McKenna W. Stanford
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
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

plt.rc('text', usetex=False)
states = cfeature.NaturalEarthFeature(category='cultural',scale='10m',facecolor='none',name='admin_1_states_provinces_shp')

file = '/glade/scratch/mckenna/AMIE/baseline/wrfout_d01_2011-12-07_12:00:00'
file2 = '/glade/scratch/mckenna/AMIE/baseline_1km/wrfout_d02_2011-12-08_00:00:00'

ncfile = xarray.open_dataset(file)
lat = np.array(ncfile['XLAT'].copy())
lon = np.array(ncfile['XLONG'].copy())
ncfile.close()

lat = np.squeeze(lat)
lon = np.squeeze(lon)

ncfile = xarray.open_dataset(file2)
lat2 = np.array(ncfile['XLAT'].copy())
lon2 = np.array(ncfile['XLONG'].copy())
ncfile.close()

lat2 = np.squeeze(lat2)
lon2 = np.squeeze(lon2)

#---------------------------
# Plot
#---------------------------
fig = plt.figure(figsize=(8.3,8.3))
Fontsize=28
ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
#ax = fig.add_subplot(111)#,projection=ccrs.PlateCarree())
states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(cfeature.OCEAN)
ax.coastlines('10m',edgecolor='black')
ax.add_feature(cfeature.LAND,edgecolor='black',color='grey')
ax.set_xlim(np.min(lon)-0.8,np.max(lon)+0.8)
ax.set_ylim(np.min(lat)-0.8,np.max(lat)+0.8)

for ii in range(400):
    ax.plot(lon[ii,0],lat[ii,0],'o',markersize=2,color='black')
    ax.plot(lon[ii,-1],lat[ii,-1],'o',markersize=2,color='black')
    ax.plot(lon[0,ii],lat[0,ii],'o',markersize=2,color='black')
    ax.plot(lon[-1,ii],lat[-1,ii],'o',markersize=2,color='black')
    
for ii in range(600):
    ax.plot(lon2[ii,0],lat2[ii,0],'o',markersize=2,color='navy')
    ax.plot(lon2[ii,-1],lat2[ii,-1],'o',markersize=2,color='navy')
    ax.plot(lon2[0,ii],lat2[0,ii],'o',markersize=2,color='navy')
    ax.plot(lon2[-1,ii],lat2[-1,ii],'o',markersize=2,color='navy')
    
x0 = np.arange(-299.5,300.5,1)
y0 = np.arange(-299.5,300.5,1)
X0,Y0 = np.meshgrid(x0,y0)
radial_dist = np.sqrt(X0**2. + Y0**2.)
#plt.contourf(lon2,lat2,radial_dist)

ax.contour(lon2,lat2,radial_dist,linewidths=[2],colors=['black'],levels=[149,151])
ax.plot(lon2[300,300],lat2[300,300],'x',color='red',markersize=12,markeredgewidth=5)

    
g1 = ax.gridlines(color='black',linestyle='dotted')
g1.xlabels_bottom = True
g1.ylabels_left = True
g1.xlabel_style = {'size':Fontsize*0.75}
g1.ylabel_style = {'size':Fontsize*0.75}
ax.set_xlabel('Longitude [$^{\\circ}$]',fontsize=Fontsize)
ax.set_ylabel('Latitude [$^{\\circ}$]',fontsize=Fontsize)

ax.text(0.1,0.85,'3 km',fontsize=Fontsize*1.2,transform=ax.transAxes,fontweight='bold')
ax.text(0.31,0.63,'1 km',fontsize=Fontsize*1.2,transform=ax.transAxes,color='navy',fontweight='bold')
ax.text(0.43,0.535,'SPOL',fontsize=Fontsize*0.75,transform=ax.transAxes,fontweight='bold')

ax.text(0.3,-0.15,'Longitude [$^{\\circ}$]',fontsize=Fontsize,transform=ax.transAxes)
ax.text(-0.15,0.35,'Latitude [$^{\\circ}$]',fontsize=Fontsize,transform=ax.transAxes,rotation=90,ha='center',va='baseline')

save_path = '/glade/u/home/mckenna/figures/amie_paper/'
outfile = save_path+'fig_01.png'
#outfile = save_path+'fig_01.eps'
plt.show()
#plt.savefig(outfile,dpi=200,bbox_inches='tight')
plt.close()
print('done')








