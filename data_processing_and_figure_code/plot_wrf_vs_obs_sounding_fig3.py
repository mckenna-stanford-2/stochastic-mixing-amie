#==========================================================================
# Title: plot_wrf_vs_obs_sounding_fig3.py
# Author: McKenna W. Stanford
# Utility: Plots a 300x300km domain-mean sounding from the WRF simulations
# and compares them to the Gan Sounings.
# Makes Fig. 3 of the manuscript.
#==========================================================================

#------------------------------------------
# Imports
#------------------------------------------
import numpy as np
import metpy
import matplotlib.pyplot as plt
import xarray
import metpy.calc as mpcalc
from metpy.plots import SkewT
from metpy.units import units
import glob
import datetime
import pickle
from netCDF4 import Dataset
import wrf

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
plt.rc('text',usetex='True')

# read in observed sounding
path1 = '/glade/scratch/mckenna/AMIE/amie_post/soundings/'
pkl_file = open(path1+'AMIE_Gan_Sounding_2011_12_07-08.p','rb')
gan_sound = pickle.load(pkl_file)
pkl_file.close()


#-------------------------------------
#-------------------------------------
# index determines the time stamp used
# for the observed soundings
#-------------------------------------
#-------------------------------------
index1 = 0 # 6-2335 UTC
index2 = 1 # 7-0235 UTC
index3 = 2 # 7-0630 UTC
index4 = 3 # 7-0835 UTC
index5 = 4 # 7-1135 UTC
index6 = 5 # 7-1435 UTC
index7 = 6 # 7-1742 UTC
index8 = 7 # 7-2035 UTC
index9 = 8 # 7-2336 UTC
index10 = 9 # 7-2336 UTC
index11 = 10 # 8-0019 UTC
index12 = 11 # 8-0235 UTC
index13 = 12 # 8-0538 UTC
index14 = 13 # 8-0620 UTC
index15 = 14 # 8-0835 UTC
index16 = 15 # 8-1137 UTC
#-------------------------------------
#-------------------------------------
#-------------------------------------
indices = [index1,index5,index10,index16]


# now read in first wrfout file (i.e. initial conditions)
path2 = '/glade/scratch/mckenna/AMIE/baseline/'

file1 = path2+'wrfout_d01_2011-12-07_00:00:00'
file2 = path2+'wrfout_d01_2011-12-07_12:00:00'
file3 = path2+'wrfout_d01_2011-12-08_00:00:00'
file4 = path2+'wrfout_d01_2011-12-08_12:00:00'

files = [file1,file2,file3,file4]

str1_1 =str(datetime.datetime(2011,12,7,0,0))
str1_2 =str(datetime.datetime(2011,12,7,12,0))
str1_3 =str(datetime.datetime(2011,12,8,0,0))
str1_4 =str(datetime.datetime(2011,12,8,12,0))

tmp_1 = datetime.datetime(2011,12,7,0,0)
tmp_2 = datetime.datetime(2011,12,7,12,0)
tmp_3 = datetime.datetime(2011,12,8,0,0)
tmp_4 = datetime.datetime(2011,12,8,12,0)



str1s = [str1_1,str1_2,str1_3,str1_4]
tmps = [tmp_1,tmp_2,tmp_3,tmp_4]

nk = len(tmps)

#---------------------------
# Start Figure
#---------------------------
fig = plt.figure(figsize=(12,12))
Fontsize=14
titles1 = ['7 Dec., 0000 UTC',\
          '7 Dec., 1200 UTC',\
          '8 Dec., 0000 UTC',\
          '8 Dec., 1200 UTC',\
         ]
titles2 = ['6 Dec., 2335 UTC',\
          '7 Dec., 1135 UTC',\
          '8 Dec., 0019 UTC',\
          '8 Dec., 1137 UTC',\
         ]

labs = ['(a)','(b)','(c)','(d)']

titdum = 0
for kk in range(nk):
    print('Sounding {}/{}'.format(kk+1,nk))
    
    
    skew = SkewT(fig,rotation=45,subplot=([2,2,kk+1]))
    
    str1 = str1s[kk]
    tmpfile = files[kk]
    index = indices[kk]
    
    ncfile = Dataset(tmpfile)
    z = wrf.getvar(ncfile,'z',meta=False).data.T
    pres = wrf.getvar(ncfile,'pres',units='hPa',meta=False).data.T
    temp = wrf.getvar(ncfile,'tc',meta=False).data.T
    u = wrf.getvar(ncfile,'ua',meta=False).data.T
    qv = wrf.getvar(ncfile,'QVAPOR',meta=False).data.T
    v = wrf.getvar(ncfile,'va',meta=False).data.T
    dp = wrf.getvar(ncfile,'td',meta=False).T
    slp = wrf.getvar(ncfile,'slp',meta=False).T
    temp2 = wrf.getvar(ncfile,'T2',meta=False).data.T
    lat = wrf.getvar(ncfile,'XLAT',meta=False).data.T
    lon = wrf.getvar(ncfile,'XLONG',meta=False).data.T
    temp2 = temp2-273.15
    dp2 = wrf.getvar(ncfile,'td2',meta=False).data.T
    #u10 = wrf.getvar(ncfile,'ua10',meta=False.data.T)
    u10_v10 = wrf.getvar(ncfile,'uvmet10',meta=False).data.T
    u10 = u10_v10[:,:,0]
    v10 = u10_v10[:,:,1]
    ncfile.close()
    
    var_dict = {'z':z,\
                'pres':pres,\
                'temp':temp,\
                'u':u,\
                'v':v,\
                'qv':qv,\
                'dp':dp,\
                'slp':slp,\
                'temp2':temp2,\
                'lat':lat,\
                'lon':lon,\
                'dp2':dp2,\
                'u10':u10,\
                'v10':v10,\
               }
    #for key,val in var_dict.items():
    #    print(key,np.shape(val),np.max(val),np.min(val))
    
    
    # To calculate Dew Point differently (?)
    idiff = 0.
    if idiff == 1.:
        ncfile = Dataset(tmpfile)
        pres3 = wrf.getvar(ncfile,'pres',units='hPa',meta=False).data.T
        qv3= wrf.getvar(ncfile,'QVAPOR',meta=False).data.T
        #dp3 = wrf.td(pres3, qv3, meta=False, units='degC')
        tdc = qv3*pres3/(0.622 + qv3)
        dp3 = (243.5*np.log(tdc) - 440.8)/(19.40 - np.log(tdc))
        dp = dp3
        ncfile.close()

    

    midpoint = 200
    radius = 50
    startid = int(midpoint-radius)
    endid = int(midpoint+radius)
    # limit to 300 km x 300 km (SPOL domain)
    z = z[startid:endid,startid:endid,:]
    pres = pres[startid:endid,startid:endid,:]
    temp = temp[startid:endid,startid:endid,:]
    dp = dp[startid:endid,startid:endid,:]
    slp = slp[startid:endid,startid:endid]
    temp2 = temp2[startid:endid,startid:endid]
    dp2 = dp2[startid:endid,startid:endid]
    u10 = u10[startid:endid,startid:endid]
    v10 = v10[startid:endid,startid:endid]


    tmp_mean_z = np.nanmean(z,axis=(0,1))
    tmp_mean_u = np.nanmean(u,axis=(0,1))
    tmp_mean_v = np.nanmean(v,axis=(0,1))
    tmp_mean_temp = np.nanmean(temp,axis=(0,1))
    tmp_mean_dp = np.nanmean(dp,axis=(0,1))
    tmp_mean_pres = np.nanmean(pres,axis=(0,1))
    tmp_mean_dp2 = np.nanmean(dp2)
    tmp_mean_temp2 = np.nanmean(temp2)
    tmp_mean_slp = np.nanmean(slp)
    tmp_mean_u10 = np.nanmean(u10)
    tmp_mean_v10 = np.nanmean(v10)
    tmp_mean_pres = np.array(tmp_mean_pres)
    tmp_mean_u10 = np.array(tmp_mean_u10)
    tmp_mean_v10 = np.array(tmp_mean_v10)
    tmp_mean_u = np.array(tmp_mean_u)
    tmp_mean_v = np.array(tmp_mean_v)

    mean_temp = np.insert(tmp_mean_temp,0,tmp_mean_temp2) * units.degC
    mean_dp = np.insert(tmp_mean_dp,0,tmp_mean_dp2) * units.degC
    mean_pres = np.insert(tmp_mean_pres,0,tmp_mean_slp) * units.hPa
    mean_z = np.insert(tmp_mean_z,0,0)
    mean_u = np.insert(tmp_mean_u,0,tmp_mean_u10)
    mean_v = np.insert(tmp_mean_v,0,tmp_mean_v10)
    mean_u = mean_u*1.94 * units.knots
    mean_v = mean_v*1.94 * units.knots
    
    #====================================
    # Plots
    #====================================
    alt = gan_sound['alt']
    mean_alt = np.nanmean(alt,axis=0)    
    wind_alts = [500,1000,1500,2000,2500,3000,3500,4000,4500,5000,\
                5500,6000,6500,7000,7500,8000,9500,10000,11000,12000,
                13000,14000,15000,16000]
    
    mean_alt = mean_alt[0:-3]
    

    
    obs_wind_alt_ids = []
    
    for wind_alt in wind_alts:
        tmpz,tmpid = find_nearest(mean_alt,wind_alt)
        obs_wind_alt_ids.append(tmpid)
        
    wrf_wind_alt_ids = []
    
    for wind_alt in wind_alts:
        tmpz,tmpid = find_nearest(mean_z,wind_alt)
        wrf_wind_alt_ids.append(tmpid)
                


    #-------------------------------------
    # Skew-T-LogP
    #-------------------------------------


    tmp_pres = gan_sound['pres'][index,:] * units.hPa
    tmp_temp = gan_sound['temp'][index,:] * units.degC
    tmp_dp = gan_sound['dp'][index,:] * units.degC
    tmp_wspd = gan_sound['wspd'][index,:] * units.meter / units.second
    tmp_wdir = gan_sound['wdir'][index,:] * units.degrees
    u, v = mpcalc.wind_components(tmp_wspd, tmp_wdir)
    u = u*1.94 * units.knots
    v = v*1.94 * units.knots
    ids = obs_wind_alt_ids
    ids2 = wrf_wind_alt_ids
    #ids = np.arange(0,202,8)
    #ids2 = np.arange(0,51,2)
    #ids2 = ids2[0:-4]




    # observed
    skew.plot(tmp_pres,tmp_temp, 'r',lw=2,label='Observed $T$')
    skew.plot(tmp_pres,tmp_dp, 'g',lw=2,label='Observed $T_{d}$')
    skew.plot(mean_pres,mean_temp, 'blue',lw=2,label = 'WRF $T$')
    skew.plot(mean_pres,mean_dp, 'orange',lw=2,label='WRF $T_{d}$')
    ##skew.plot_barbs(mean_pres,mean_u,mean_v)
    skew.plot_barbs(tmp_pres[ids],u[ids],v[ids],xloc=1.)
    #tmppres = np.ma.getdata(mean_pres)
    #tmpu = np.ma.getdata(mean_u)
    #tmpv = np.ma.getdata(mean_v)
    
    
    #skew.plot_barbs(mean_pres[ids2],mean_u[ids2],mean_v[ids2],xloc=1.12,color='blue')
    skew.plot_barbs(mean_pres[ids2],mean_u[ids2],mean_v[ids2],xloc=1.09,color='blue')

    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-40, 40)

    skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)

    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

    plt.xlabel('Temperature [$^{o}$C]',fontsize=Fontsize)
    plt.ylabel('Pressure [hPa]',fontsize=Fontsize)
    plt.tick_params(labelsize=Fontsize)
    str2 = str(gan_sound['time'][index])

    tmpyear = str(tmps[kk].year)
    tmpmonth = str(tmps[kk].month)
    tmpday = str(tmps[kk].day)
    tmphour = str(tmps[kk].hour)
    tmpminute = str(tmps[kk].minute)

    heights = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]) * units.km
    std_pressures = mpcalc.height_to_pressure_std(heights)
    for height_tick, p_tick in zip(heights, std_pressures):
        trans, _, _ = skew.ax.get_yaxis_text1_transform(0)
        skew.ax.text(0.02, p_tick, '---{:~d}'.format(height_tick), transform=trans,fontsize=Fontsize*0.7)
        
    trans, _, _ = skew.ax.get_yaxis_text1_transform(0)    
    skew.ax.text(-0.3,80,labs[titdum],fontweight='bold',transform=trans,fontsize=Fontsize*3)

        
    plt.title('WRF Average: '+titles1[titdum]+'\nObserved Gan Island: '+titles2[titdum],\
              fontsize=Fontsize*1.5,va='bottom')
    
    titdum+=1
        
plt.legend(fontsize=Fontsize*1.5,bbox_to_anchor=(0.5,-0.1),ncol=2,framealpha=False)

plt.subplots_adjust(wspace=0.425,hspace=0.375)

savepath = '/glade/u/home/mckenna/figures/amie_paper/'
outfile = 'fig_03.png'
#outfile = 'fig_03.eps'
plt.savefig(savepath+outfile,dpi=260,bbox_inches='tight')
#plt.show()
plt.close()


