#==========================================================================
# Title: compute_wrf_inst_rain_rates_1km.py
# Author: McKenna W. Stanford
# Utility: Reads in dictionaries named baseline_1km_rain_vars_at_2.5-km.p
# and computes the instantaneous rain rates for 12 UTC on 7 Dec. through
# 12 UTC on 8 Dec.
# Writes this to a dictionary named '<SIM_NAME>_inst_rain_rates.p.
#==========================================================================

import numpy as np
import pickle
import xarray
import scipy.special as scs
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature
import matplotlib
import glob
from netCDF4 import Dataset
from matplotlib.cm import get_cmap


#------------------------------------------------------------
# Function that computes and returns rain size distribution 
# parameters. Pulled from P3.
#-----------------------------------------------------------
def get_rain_dsd(tmpqr,tmpnr):

    # in : qr,nr
    # inout: nr_grid (in return statement)
    # out: lamr,cdistr,logn0r (in return statement)

    if tmpqr >= qsmall:
        # find spot in lookup table
        # (scaled N/q for lookup table in parameter space)
        tmpnr = max(tmpnr,nsmall)
        inv_dum = (tmpqr/(cons1*tmpnr*6.))**thrd

        lamr = (cons1*tmpnr*(mu_r+3.)*(mu_r+2.)*(mu_r+1.)/(tmpqr))**thrd # recalculate slope based on mu_r
        lammax = (mu_r+1.)*1.e+5 # check for slope
        lammin = (mu_r+1.)*1250. # set to small value since breakup is explicitly included (mean size 0.8 mm)

        # apply lambda limiters for rain
        if (lamr < lammin):
            lamr = lammin
            tmpnr = np.exp(3.*np.log(lamr)+np.log(tmpqr)+np.log(scs.gamma(mu_r+1.))-np.log(scs.gamma(mu_r+4.)))/(cons1)
        elif (lamr > lammax):
            lamr = lammax
            tmpnr = np.exp(3.*np.log(lamr)+np.log(tmpqr)+np.log(scs.gamma(mu_r+1.))-np.log(scs.gamma(mu_r+4.)))/(cons1)


        logn0r = np.log10(tmpnr)+(mu_r+1.)*np.log10(lamr)-np.log10(scs.gamma(mu_r+1.)) # note: logn0r is calcualted as log10(n0r)
        cdistr = tmpnr/scs.gamma(mu_r+1.)

    else:
        lamr = 0.
        cdistr = 0.
        logn0r = 0.

    return tmpnr,lamr,cdistr,logn0r



#------------------------------------------------------------
# Function that finds indices in a lookup table. Pulled from P3.
#-----------------------------------------------------------
def find_lookupTable_indices_3(mu_r,tmplamr):
    # in: mu_r, lamr
    # out: dumii,dumjj, dum1, rdumii, rdumjj, inv_dum3

    # find location in scaled mean size space
    dum1 = (mu_r+1.)/tmplamr
    
    if (dum1 <= 195.e-6):
        inv_dum3 = 0.1
        
        rdumii = (dum1*1.e6+5.)*inv_dum3
        rdumii = rdumii-1 # for python
        rdumii = max(rdumii,0.)
        rdumii = min(rdumii,19.)
        
        dumii = int(rdumii)
        dumii = max(dumii,0)
        dumii = min(dumii,19)
        
    elif (dum1 > 195.e-6):
        
        inv_dum3 = thrd*0.1
        rdumii = (dum1*1.e+6-195.)*inv_dum3 + 20.
        rdumii = rdumii - 1 #for python
        rdumii = max(rdumii, 19.)
        rdumii = min(rdumii,299.)
        
        dumii = int(rdumii)
        dumii = max(dumii,19)
        dumii = min(dumii,298)
    
    
    #find location in mu_r space
    rdumjj = mu_r+1.
    rdumjj = rdumjj-1 # for python
    rdumjj = max(rdumjj,0.)
    rdumjj = min(rdumjj,9.)
    dumjj = int(rdumjj)
    dumjj = max(dumjj,0)
    dumjj = min(dumjj,8)

    return dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3	

#------------------------------------------------------------
# Define local variables
#------------------------------------------------------------
mu_r = 0.
rd = 287.15
qsmall = 1.e-14
pi = 3.14159265
rhow = 997.
sxth = 1./6.
thrd = 1./3.
piov6 = pi*sxth
nsmall = 1.e-16
cons1 = piov6*rhow
rhosur = 100000./(rd*273.15)

# Read in lookup table
pkl_file = open('vm_table','rb')
vm_table = pickle.load(pkl_file)
pkl_file.close()


path = '/glade/scratch/mckenna/AMIE/amie_post/'
pkl_file = open(path+'rain_vars/baseline_1km_rain_vars_at_2.5-km.p','rb')
rain_vars_dict = pickle.load(pkl_file)
pkl_file.close()


#if True:
if False:
    for key,val in rain_vars_dict.items():
        print(key,np.shape(val))


rhofacr = (rhosur/rain_vars_dict['rho_air'][24:])**0.54   

qr = rain_vars_dict['qr'][24:,:,:]
nr = rain_vars_dict['nr'][24:,:,:]
rho_air = rain_vars_dict['rho_air'][24:,:,:]
time = rain_vars_dict['time_arr'][24:]

nt = len(qr[:,0,0])
nx = len(qr[0,:,0])
ny = len(qr[0,0,:])


rain_rate = np.zeros((nt,nx,ny))

#for tt in range(-2,-1):
for tt in range(nt):
    print('Time step {}/{}:'.format(tt+1,nt))
    for ii in range(nx):
        for jj in range(ny):

            if qr[tt,ii,jj] <= qsmall:
                rain_rate[tt,ii,jj] = 0.0
                continue
            else:
                pass

            # Grab DSD parameters
            data = get_rain_dsd(qr[tt,ii,jj],nr[tt,ii,jj])
            tmp_nr = data[0]
            tmp_lamr = data[1]
            tmp_cdistr= data[2]
            tmp_logn0r= data[3]

            data2 = find_lookupTable_indices_3(mu_r,tmp_lamr)
            dumii = data2[0]
            dumjj = data2[1]
            dum1 = data2[2]
            rdumii = data2[3]
            rdumjj = data2[4]
            inv_dum3 = data2[5]

            # mass-weighted fall speed
            dum1 = vm_table[dumii,dumjj]+(rdumii-float(dumii))\
                *(vm_table[dumii+1,dumjj]-vm_table[dumii,dumjj])

            dum2 = vm_table[dumii,dumjj+1]+(rdumii-float(dumii))\
                *(vm_table[dumii+1,dumjj+1]-vm_table[dumii,dumjj+1])

            #interpolated:
            Vt_qr = dum1 + (rdumjj-float(dumjj))*(dum2-dum1)
            Vt_qr = Vt_qr*rhofacr[tt,ii,jj]
            rain_rate[tt,ii,jj] = Vt_qr*qr[tt,ii,jj]*1000.*3600.*rho_air[tt,ii,jj]/997.

#--------------------------------------------------
# Plot rain rates at 2.5 km AGL (optional)
#--------------------------------------------------
iplot = False

if iplot:
    lat = rain_vars_dict['lat']
    lon = rain_vars_dict['lon']
    time_id = -2
    fig = plt.figure(figsize=(8,8))
    Fontsize=14
    ax = fig.add_subplot(111,projection=crs.PlateCarree())
    levs = [0.01,0.05,0.1,0.5,1,5,10,50,100]
    tmpplot = ax.contourf(lon,lat,rain_rate[time_id,:,:],levels=levs,\
                        cmap=get_cmap('nipy_spectral'),\
                        transform=crs.PlateCarree(),\
                        norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256))

    cbar_ax = fig.add_axes([0.13,0.0,0.75,0.04])   
    # Set color mappable
    range_min = 0
    range_max = 60
    cmap = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.BoundaryNorm(boundaries=levs,ncolors=256),\
        cmap = plt.get_cmap('nipy_spectral'))    
    cmap.set_array([])
    cbar = fig.colorbar(cmap,orientation='horizontal',cax=cbar_ax)
    cbar.ax.tick_params(labelsize=Fontsize)
    cbar.set_label('Instantaneous Rain Rate [mm hr$^{-1}$]',fontsize=Fontsize)
    ax.tick_params(labelsize=Fontsize)
    g1 = ax.gridlines(color='black',linestyle='solid',lw=1.)
    g1.xlabels_bottom = True
    g1.ylabels_left = True
    g1.xlabel_style = {'size':Fontsize}
    g1.ylabel_style = {'size':Fontsize}
    ax.set_aspect('auto')
    ax.set_xlim(np.min(lon),np.max(lon))
    ax.set_ylim(np.min(lat),np.max(lat))
    time_str = time[time_id].strftime('%Y/%m/%d %H:%M UTC')
    ax.set_title('2.5 km AGL, {}'.format(time_str),fontsize=Fontsize*1.5)

    plt.show()
    plt.close()

    #-------------------------
    # End Plotting
    #-------------------------

out_dict = {'rain_rate':rain_rate,\
              'lon':rain_vars_dict['lon'],\
              'lat':rain_vars_dict['lat'],\
              'time':time}
# Print diagnostics
#if False:
if True:
    for key,val in out_dict.items():
        print(key,np.shape(val),np.min(val),np.max(val))

print('Writing to dictionary...')
savdir = path+'inst_rain_rates/'
f = open(savdir+'baseline_1km_inst_rain_rates.p','wb')
pickle.dump(out_dict,f)
f.close()    
    
print('Completed writing to dictinary. Done.')



    
