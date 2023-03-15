#==========================================================================
# Title: compress_amie_gan_soundings.py
# Author: McKenna W. Stanford
# Utility: Compresses soundings from AMIE Gan Island to a single file
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


path = '/glade/scratch/mckenna/AMIE/soundings/'

ncfile = xarray.open_dataset(path+'DYNAMO_43599_L4_5hPa_Sounding_Data.nc')
ncfile.close()

n_soundings = np.array(ncfile['n_soundings'].copy())
date = np.array(ncfile['release_date_enc'])
time_utc = np.array(ncfile['release_time_enc'])
pres = np.array(ncfile['p'].copy()) #hPa
temp = np.array(ncfile['T'].copy()) #degC
dp = np.array(ncfile['Td'].copy()) #degC
wspd = np.array(ncfile['wind_spd'].copy()) #m/s
wdir = np.array(ncfile['wind_dir'].copy()) #degree
lon = np.array(ncfile['lon'].copy()) #deg. E
lat = np.array(ncfile['lat'].copy()) #deg. N
alt = np.array(ncfile['alt'].copy()) #alt. above MSL

nt = len(pres[:,0])


#==================================================
# Make datetime array
#==================================================
year = []
month = []
day = []
hour = []
minute = []
second = []

time = []
for tt in range(nt):
    tmp_date = str(date[tt])
    
    tmp_year_str = int(tmp_date[0:4])
    year.append(tmp_year_str)
    tmp_month_str = int(tmp_date[4:6])
    month.append(tmp_month_str)
    tmp_day_str = int(tmp_date[6:8])
    day.append(tmp_day_str)
    
    tmp_utc_time = str(time_utc[tt])
    if len(tmp_utc_time) == 3:
        hour.append(0)
        minute.append(int(tmp_utc_time[0]))
        second.append(int(tmp_utc_time[1:2]))
    elif len(tmp_utc_time) == 4:   
        hour.append(0)
        minute.append(int(tmp_utc_time[0:2]))
        second.append(int(tmp_utc_time[2:4]))
    elif len(tmp_utc_time) == 5:   
        hour.append(int(tmp_utc_time[0:1]))
        minute.append(int(tmp_utc_time[1:3]))
        second.append(int(tmp_utc_time[3:5]))
    else:
        hour.append(int(tmp_utc_time[0:2]))
        minute.append(int(tmp_utc_time[2:4]))
        second.append(int(tmp_utc_time[4:6]))        
        
    time.append(datetime.datetime(year[tt],month[tt],day[tt],hour[tt],minute[tt]))
    
# limit to only cover 7th and 8th december    
time = time[659:678]
time_utc = time_utc[659:678]
pres = pres[659:678,:]
temp = temp[659:678,:]
dp = dp[659:678,:]
alt = alt[659:678,:]
lon = lon[659:678,:]
lat = lat[659:678,:]
wspd = wspd[659:678,:]
wdir = wdir[659:678,:]
date = date[659:678]

# save these to a dictionary
obs_dict = {'time':time,'time_utc':time_utc,'pres':pres,'temp':temp,'dp':dp,\
           'alt':alt,'lon':lon,'lat':lat,'wspd':wspd,'wdir':wdir,'date':date}


savdir = '/glade/scratch/mckenna/AMIE/amie_post/soundings/'
f = open(savdir+'AMIE_Gan_Sounding_2011_12_07-08.p','wb')
pickle.dump(obs_dict,f)
f.close()





