#==========================================================================
# Title: degrade_1km_to_3km.py
# Author: McKenna W. Stanford
# Utility: Coarse-grains 1-km simulation to 3-km grid spacing for all
# relevant variables needed to perform analyses. These are saved for
# individual times and written to NetCDF files.
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
import pandas as pd
import datetime as dt


path = '/glade/scratch/mckenna/AMIE/baseline_1km/'
files = sorted(glob.glob(path+'wrfout_d02*'))
# Limit to after 00 UTC
files = files[24:]

#files = files[0:25]
#files = files[25:50]
#files = files[50:75]
files = files[75:]
nt = len(files)

for tt in range(nt):
    print('Time step {}/{}:'.format(tt+1,nt))  
    print('Reading in file:',files[tt])
    ncfile = Dataset(files[tt])
    lat = wrf.getvar(ncfile,'lat',meta=False).data.T
    lon = wrf.getvar(ncfile,'lon',meta=False).data.T    
    z = wrf.getvar(ncfile,'height',units='m',meta=False).data.T    
    pressure = wrf.getvar(ncfile,'pres',units='hPa',meta=False).data.T    
    w = wrf.getvar(ncfile,'wa',meta=False).data.T    
    u = wrf.getvar(ncfile,'ua',meta=False).data.T    
    v = wrf.getvar(ncfile,'va',meta=False).data.T    
    qice = wrf.getvar(ncfile,'QICE',meta=False).data.T    
    qrain = wrf.getvar(ncfile,'QRAIN',meta=False).data.T    
    qcloud = wrf.getvar(ncfile,'QCLOUD',meta=False).data.T    
    tracer = wrf.getvar(ncfile,'TR1',meta=False).data.T    
    hdi = wrf.getvar(ncfile,'H_DIABATIC',meta=False).data.T    
    cconevp = wrf.getvar(ncfile,'CCONEVP',meta=False).data.T    
    rconevp = wrf.getvar(ncfile,'RCONEVP',meta=False).data.T    
    qir = wrf.getvar(ncfile,'QIR',meta=False).data.T    
    qib = wrf.getvar(ncfile,'QIB',meta=False).data.T    
    dbz= wrf.getvar(ncfile,'REFL_10CM',meta=False).data.T      
    tv= wrf.getvar(ncfile,'tv',units='K',meta=False).data.T  
    time = pd.to_datetime(wrf.getvar(ncfile,'times').data)
    ncfile.close()
    print('Completed reading in file...')
    
    nx = len(pressure[:,0,0])
    ny = len(pressure[0,:,0])
    nz = len(pressure[0,0,:])
    
    ze = np.zeros((nx,ny,nz))
    tmpid = np.where(dbz > -99.)
    ze[tmpid] = 10.**(dbz[tmpid]/10.)
    tmpid = np.where(dbz == -99)
    ze[tmpid] = 0.
    
    rd = 287.04
    rho_air = 100.*pressure/(rd*tv) #kg/m^3
    
    nx = len(tv[:,0,0])
    ny = len(tv[0,:,0])
    nz = len(tv[0,0,:])
    
    fac = 3.
    nx_3km = int(np.floor(nx/fac))
    ny_3km = int(np.floor(ny/fac))
    
    tracer_3km = np.zeros((nx_3km, ny_3km,nz))
    pressure_3km = np.zeros((nx_3km, ny_3km,nz))
    qice_3km = np.zeros((nx_3km, ny_3km,nz))
    qrain_3km = np.zeros((nx_3km, ny_3km,nz))
    qcloud_3km = np.zeros((nx_3km, ny_3km,nz))
    qir_3km = np.zeros((nx_3km, ny_3km,nz))
    qib_3km = np.zeros((nx_3km, ny_3km,nz))
    w_3km = np.zeros((nx_3km, ny_3km,nz))
    u_3km = np.zeros((nx_3km, ny_3km,nz))
    v_3km = np.zeros((nx_3km, ny_3km,nz))
    hdi_3km = np.zeros((nx_3km, ny_3km,nz))
    cconevp_3km = np.zeros((nx_3km, ny_3km,nz))
    rconevp_3km = np.zeros((nx_3km, ny_3km,nz))
    z_3km = np.zeros((nx_3km, ny_3km,nz))
    tv_3km = np.zeros((nx_3km, ny_3km,nz))
    ze_3km = np.zeros((nx_3km, ny_3km,nz))
    dbz_3km = np.zeros((nx_3km, ny_3km,nz))
    lon_3km = np.zeros((nx_3km, ny_3km))
    lat_3km = np.zeros((nx_3km, ny_3km))
    
    
    print('Coarse graining...')
    fac = int(fac)
    # loop through heights
    for kk in range(nz):
        print('Height {}/{}:'.format(kk+1,nz))
        for ii in range(nx_3km):
            for jj in range(ny_3km):
                if kk == 0:
                    lon_3km[ii,jj] = np.mean(lon[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac)])
                    lat_3km[ii,jj] = np.mean(lat[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac)])
                    
                w_3km[ii,jj,kk] = np.mean(w[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                u_3km[ii,jj,kk] = np.mean(u[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                v_3km[ii,jj,kk] = np.mean(v[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                hdi_3km[ii,jj,kk] = np.mean(hdi[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                cconevp_3km[ii,jj,kk] = np.mean(cconevp[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                rconevp_3km[ii,jj,kk] = np.mean(rconevp[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                tracer_3km[ii,jj,kk] = np.mean(tracer[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                ze_3km[ii,jj,kk] = np.mean(ze[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                qrain_3km[ii,jj,kk] = np.mean(qrain[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                qcloud_3km[ii,jj,kk] = np.mean(qcloud[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                qice_3km[ii,jj,kk] = np.mean(qice[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                qir_3km[ii,jj,kk] = np.mean(qir[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                qib_3km[ii,jj,kk] = np.mean(qib[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                pressure_3km[ii,jj,kk] = np.mean(pressure[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                tv_3km[ii,jj,kk] = np.mean(tv[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                z_3km[ii,jj,kk] = np.mean(z[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
                ze_3km[ii,jj,kk] = np.mean(ze[ii*fac:((ii*fac)+fac),jj*fac:((jj*fac)+fac),kk])
    
    print('Completed coarse graining.')
    tmpid = np.where(ze_3km == 0.)
    dbz_3km[tmpid] = -999.
    tmpid = np.where(ze_3km > 0.)
    dbz_3km[tmpid] = 10.*np.log10(ze_3km[tmpid])
    
    rho_air_3km = 100.*pressure_3km/(rd*tv_3km)
    qt_3km = qice_3km+qrain_3km+qcloud_3km
    twc_3km = qt_3km*rho_air_3km*1.e3

    
    
    print('Writing to files...')
    old_path = files[tt]
    tags = str.split(old_path,'/')
    old_file_name = tags[-1]
    outfile = path+'3km/'+old_file_name+'_1km_coarse_grained_to_3km.nc'
    print('Outfile name:',outfile)    
    ncfile = netcdf.netcdf_file(outfile,'w',version=2)
    ncfile.createDimension('west_east',nx_3km)
    ncfile.createDimension('north_south',ny_3km)
    ncfile.createDimension('bottom_top',nz)    
    ncfile.createDimension('time_len',1)    

    tmp_w = ncfile.createVariable('w',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_u = ncfile.createVariable('u',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_v = ncfile.createVariable('v',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_z = ncfile.createVariable('z',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_dbz = ncfile.createVariable('dbz',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_qice = ncfile.createVariable('qice',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_qrain = ncfile.createVariable('qrain',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_qcloud = ncfile.createVariable('qcloud',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_hdi = ncfile.createVariable('hdi',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_qir = ncfile.createVariable('qir',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_qib = ncfile.createVariable('qib',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_tracer = ncfile.createVariable('tracer',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_pressure = ncfile.createVariable('pressure',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_rho_air = ncfile.createVariable('rho_air',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_tv = ncfile.createVariable('tv',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_cconevp = ncfile.createVariable('cconevp',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_rconevp = ncfile.createVariable('rconevp',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_twc = ncfile.createVariable('twc',dimensions=['west_east','north_south','bottom_top'],type='float64')
    tmp_lon = ncfile.createVariable('lon',dimensions=['west_east','north_south'],type='float64')
    tmp_lat = ncfile.createVariable('lat',dimensions=['west_east','north_south'],type='float64')
    tmp_time = ncfile.createVariable('time',dimensions=['time_len',],type='int32')
    
    tmp_w[:] = w_3km
    tmp_u[:] = u_3km
    tmp_v[:] = v_3km
    tmp_z[:] = z_3km
    tmp_dbz[:] = dbz_3km
    tmp_twc[:] = twc_3km
    tmp_qice[:] = qice_3km
    tmp_qrain[:] = qrain_3km
    tmp_qcloud[:] = qcloud_3km
    tmp_qir[:] = qir_3km
    tmp_qib[:] = qib_3km
    tmp_tv[:] = tv_3km
    tmp_hdi[:] = hdi_3km
    tmp_cconevp[:] = cconevp_3km
    tmp_rconevp[:] = rconevp_3km
    tmp_tracer[:] = tracer_3km
    tmp_pressure[:] = pressure_3km
    tmp_rho_air[:] = rho_air_3km
    tmp_lon[:] = lon_3km
    tmp_lat[:] = lat_3km

    # Make epoch time stamp
    epoch_time = int(pd.Timestamp(time).timestamp())
    #print(epoch_time)
    #print(pd.to_datetime(epoch_time,unit='s'))
    tmp_time[:] = epoch_time
    
    ncfile.close()
    print('Completed writing. Done.')
    
    
    
    
    
    
 