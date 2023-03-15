#==========================================================================
# Title: plot_stochastic_patterns_fig2.py
# Author: McKenna W. Stanford
# Utility: Plots horizontal cross sections of stochastic patterns.
# Produces Fig. 2 of the manuscript.
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
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.offsetbox import AnchoredText
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

sav_path = '/glade/u/home/mckenna/figures/amie_paper/'
plt.rc('text', usetex=True)
#------------------------------------------
# Parameters
#------------------------------------------

#------------------------------------------
# Functions
#------------------------------------------
# make function that finds the nearest
# element in array to the given value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

path_short = '/glade/scratch/mckenna/AMIE/stoch1_short/'
path_long = '/glade/scratch/mckenna/AMIE/stoch1_long/'

dumid = 110
files_short = sorted(glob.glob(path_short+'/wrfout*'))[dumid:]
files_long = sorted(glob.glob(path_long+'/wrfout*'))[dumid:]

nt = len(files_short)

titles = ['STOCH1\_SHORT','STOCH1\_LONG']

#=====================================================
# Plot Fig. 1 of manuscript
#=====================================================

# Will loop through and find one we like
for tt in range(nt):
    print('Time step {}/{}:'.format(tt+1,nt))
    print('tt:',tt+60)
    #cmap = 'seismic'
    cmap = 'BrBG'
    
    # Start Figure
    Fontsize=24
    fig = plt.figure(figsize=(8.3,16))
    ax1 = fig.add_subplot(211,projection=ccrs.PlateCarree())
    ax2 = fig.add_subplot(212,projection=ccrs.PlateCarree())
    
    axlist = [ax1,ax2]
    
    levs = np.arange(-6,6.25,0.25)
    
    #STOCH SHORT
    ncfile = xarray.open_dataset(files_short[tt])
    rand_pert_short =ncfile['RAND_PERT'].values.squeeze().T
    lon =ncfile['XLONG'].values.squeeze().T
    lat =ncfile['XLAT'].values.squeeze().T
    ncfile.close()


    tmp_plot = ax1.contourf(lon,lat,rand_pert_short[:,:,0],levels=levs,\
                                    norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256),\
                                       cmap=get_cmap(cmap),\
                                   transform=ccrs.PlateCarree(),extend='both')
    #STOCH LONG
    ncfile = xarray.open_dataset(files_long[tt])
    rand_pert_long =ncfile['RAND_PERT'].values.squeeze().T
    lon =ncfile['XLONG'].values.squeeze().T
    lat =ncfile['XLAT'].values.squeeze().T
    ncfile.close()

    tmp_plot = ax2.contourf(lon,lat,rand_pert_long[:,:,0],levels=levs,\
                                    norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256),\
                                       cmap=get_cmap(cmap),\
                                   transform=ccrs.PlateCarree(),extend='both')

        
    for ax in axlist:
        ax.coastlines('10m')
        ax.set_xlim(np.min(lon),np.max(lon))
        ax.set_ylim(np.min(lat),np.max(lat))    
        g1 = ax.gridlines(color='black',linestyle='dotted',lw=2)
        g1.xlabels_bottom = True
        g1.ylabels_left = True
        g1.xlabel_style = {'size':Fontsize}
        g1.ylabel_style = {'size':Fontsize}
        ax.set_aspect('auto') 
    
        
    plt.text(0.,1.05,'(a) STOCH\_SHORT',fontsize=Fontsize*1.5,fontweight='bold',\
         transform=ax1.transAxes,ha='left')
    plt.text(0.,1.05,'(b) STOCH\_LONG',fontsize=Fontsize*1.5,fontweight='bold',\
         transform=ax2.transAxes,ha='left')
    
    plt.text(-0.15,0.5,'Latitude [$^{\\circ}$]',fontsize=Fontsize,\
         transform=ax1.transAxes,ha='center',rotation=90,va='center')
    plt.text(-0.15,0.5,'Latitude [$^{\\circ}$]',fontsize=Fontsize,\
         transform=ax2.transAxes,ha='center',rotation=90,va='center')
    plt.text(0.5,-0.13,'Longitude [$^{\\circ}$]',fontsize=Fontsize,\
         transform=ax1.transAxes,ha='center',va='center')
    plt.text(0.5,-0.13,'Longitude [$^{\\circ}$]',fontsize=Fontsize,\
         transform=ax2.transAxes,ha='center',va='center')    
    
    
    cbar_ax = fig.add_axes([0.13,0.02,0.75,0.035])   
    # Set color mappable
    range_min = -6
    range_max = 6
    cmap = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256),\
    cmap = plt.get_cmap(cmap))    
    cmap.set_array([])   
    
    dum_levs = levs[::8]
    cbar = fig.colorbar(cmap,orientation='horizontal',cax=cbar_ax,\
                       ticks=dum_levs)
    cbar.ax.tick_params(labelsize=Fontsize)
    cbar.set_label('Perturbation Magnitude',fontsize=Fontsize)        
        
    plt.subplots_adjust(hspace=0.325)


    outfile = 'fig_02.png'
    plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
    #plt.show()
    plt.close()
    print(aaa)

