# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 17:29:27 2021

@author: tilloal
"""

from nco import Nco
nco = Nco()
import re
import xarray as xr
import numpy as np
from netCDF4 import Dataset
import time
import pandas as pd
from xarray import concat
# Import the os module
import os

#%%

## Script 2/4

# script that creates yearly files of daily ERA5* variables and from aggregated hourly values
# and creates yealry files from monthly files in the specified years.
 
# Requires:
# 1) monthly files of hourly data for the specified variables in the specified years
# (Outputs from "01_ERA5land_downloader.py")
 
# Output:
# 1) separate netCDF file for chosen daily variable for each year

#%%
# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))
#%%
# key text file containing (1) proxy key, (2) CDS API key, (3) working directory of the data
keys=[]
with open('Keys.txt') as f:
    for line in f:            
        keys.append(line) 
        
# Change the current working directory
os.chdir(keys[2])
#%%
# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))


# Parameters for the generation of new netcdf files
__OUTPUT_FILE_EXT = '.nc'

__NETCDF_DATASET_FORMAT = 'NETCDF4_CLASSIC'
__NETCDF_CONVENTIONS = 'CF-1.6'
__NETCDF_SOURCE_SOFTWARE = 'Python netCDF4'

__NETCDF_VAR_TIME_DIMENSION = None
__NETCDF_VAR_TIME_CALENDAR_TYPE = 'proleptic_gregorian'

__NETCDF_VAR_DATA_TYPE = 'f8'
__NETCDF_VALUE_DATA_TYPE = 'f4'
__NETCDF_COORDINATES_DATA_TYPE = 'i4'

__KEY_STANDARD_NAME = 'value_standard_name'
__KEY_LONG_NAME = 'value_long_name'
__KEY_UNIT = 'value_unit'
__KEY_OFFSET = 0
__KEY_SCALE_FACTOR = 1
__KEY_VMIN = -400
__KEY_VMAX = 400

# compressions factors for variables tp, ws, rg, rn, td
factorz=pd.read_csv('compression_factors.csv')

#%%
# function to estimate optimal add_offset and scale_factor for other variables
def compute_scale_and_offset(min, max, n):
    # stretch/compress data to the available packed range
    scale_factor = (max - min) / (2 ** n - 1)
    # translate the range to be symmetric about zero
    add_offset = min + 2 ** (n - 1) * scale_factor
    return (scale_factor, add_offset)



__meteo_vars_config = {
    'tn' : {__KEY_UNIT : 'celcius', __KEY_STANDARD_NAME : 'tn', __KEY_LONG_NAME : 'min_temperature',
            __KEY_OFFSET : 0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : -700, __KEY_VMAX : 700},
    'tp' : {__KEY_UNIT : 'mm', __KEY_STANDARD_NAME : 'tp', __KEY_LONG_NAME : 'total_precipitation',
            __KEY_OFFSET : factorz['add_offset'].values[5], __KEY_SCALE_FACTOR : factorz['scale_factor'].values[5], __KEY_VMIN : -1, __KEY_VMAX : 7000},
    'ta' : {__KEY_UNIT : 'celcius', __KEY_STANDARD_NAME : 'ta', __KEY_LONG_NAME : 'mean_temperature',
            __KEY_OFFSET : factorz['add_offset'].values[1], __KEY_SCALE_FACTOR : factorz['scale_factor'].values[1], __KEY_VMIN : -700, __KEY_VMAX : 700},
    'td' : {__KEY_UNIT : 'celcius', __KEY_STANDARD_NAME : 'td', __KEY_LONG_NAME : 'mean_dewpoint_temperature',
            __KEY_OFFSET : factorz['add_offset'].values[4], __KEY_SCALE_FACTOR : factorz['scale_factor'].values[4], __KEY_VMIN : -700, __KEY_VMAX : 700},
    'tx' : {__KEY_UNIT : 'celcius', __KEY_STANDARD_NAME : 'tx', __KEY_LONG_NAME : 'max_temperature',
            __KEY_OFFSET : 0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : -400, __KEY_VMAX : 400},
    'ws' : {__KEY_UNIT : 'm/s', __KEY_STANDARD_NAME : 'ws', __KEY_LONG_NAME : 'avg_wind_speed',
            __KEY_OFFSET : factorz['add_offset'].values[0], __KEY_SCALE_FACTOR : factorz['scale_factor'].values[0], __KEY_VMIN : -400, __KEY_VMAX : 400},
    'rg' : {__KEY_UNIT : 'J/m2/d', __KEY_STANDARD_NAME : 'ssr', __KEY_LONG_NAME : 'surface_net_solar_radiation',
            __KEY_OFFSET : factorz['add_offset'].values[2], __KEY_SCALE_FACTOR : factorz['scale_factor'].values[2]},
    'rn' : {__KEY_UNIT : 'J/m2/d', __KEY_STANDARD_NAME : 'str', __KEY_LONG_NAME : 'surface_net_thermal_radiation',
            __KEY_OFFSET : factorz['add_offset'].values[3], __KEY_SCALE_FACTOR : factorz['scale_factor'].values[3]},
}

vars_list = ' | '.join(__meteo_vars_config.keys())
#%



#%%

# function to compute mean wind speed from u and v components of wind
def wind_uv_to_spd(U,V):
    """
    Calculates the wind speed from the u and v wind components
    Inputs:
      U = west/east direction (wind from the west is positive, from the east is negative)
      V = south/noth direction (wind from the south is positive, from the north is negative)
    """
    WSPD=np.sqrt(U**2+V**2)
    return WSPD

# function to write a new netcdf file 
def writenetcdf(yr, vr, dt, scale_factor, add_offset):   
        rfile = dt.copy()
        tziz=(rfile['time'])
        d1=(tziz[0])
        d2=np.datetime_as_string(d1, unit='D')
        d3=" 00:00:00"
        d4=d2+d3
        tunits= 'days since ' + d2
        d1=os.getcwd() + '\\'
        d2=  vr + '\\0.1_deg\e5ldx_' + vr + '_'+ yr + '.nc'
        namenc=d1+d2
        nf2=Dataset(namenc,mode='w',format='NETCDF4_CLASSIC') 
        nf2.history = 'Created Nov 2021' #####
        nf2.Conventions = 'CF-1.6'
        nf2.Source_Software = 'Python netCDF4'
        nf2.reference = 'A global daily high-resolution gridded meteorological data set for 1979-2019'  #####
        nf2.title = 'Lisflood meteo maps 1981 for EUROPE setting Nov. 2021'
        nf2.keywords = 'Lisflood, Global'
        nf2.source = 'ERA5-land'
        nf2.institution = 'European Commission - Economics of climate change Unit (JRC.C.6) : https://ec.europa.eu/jrc/en/research-topic/climate-change'
        nf2.comment = 'The timestamp marks the end of the aggregation interval for a given map.'
         
         #Dimension
        nf2.createDimension('lon', 756)
        nf2.createDimension('lat', 501)
        nf2.createDimension('time', None)
         #Variables
        longitude = nf2.createVariable('lon','f4',('lon',), complevel=4, zlib=True)  ###('lon','f8',('lon',))
        longitude.standard_name= 'Longitude'
        longitude.long_name= 'Longitude'
        longitude.units =  'degrees_east'
        
        latitude = nf2.createVariable('lat','f4',('lat',), complevel=4, zlib=True)  ###('lon','f8',('lon',))
        latitude.standard_name= 'Latitude'
        latitude.long_name= 'Latitude'
        latitude.units = 'degrees_north'
        
        time = nf2.createVariable('time', 'i4', ('time',), complevel=4, zlib=True)
        time.standard_name = 'time'
        time.units = tunits
        time.frequency = '1'
        time.calendar = 'proleptic_gregorian'
        
        proj = nf2.createVariable('wsg_1984', 'i4')
        proj.grid_mapping_name = 'latitude_longitude'
        #proj.false_easting= ''
        #proj.false_northing= ''
        #proj.longitude_of_projection_origin= ''
        #proj.latitude_of_projection_origin= ''
        proj.semi_major_axis= '6378137.0'
        proj.inverse_flattening='298.257223563'
        proj.proj4_params='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
        proj.EPSG_code='EPSG:4326'
        #proj.spatial_ref='GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]]'
        
        if vr in ['ssr','str']:
            E5landata = nf2.createVariable(vr, 'i2',('time','lat','lon',), zlib=True, complevel=4, fill_value=-32767)
        else:
            E5landata = nf2.createVariable(vr, 'i2',('time','lat','lon',), zlib=True, complevel=4, fill_value=-32767)
        
        E5landata.standard_name = __meteo_vars_config[vr][__KEY_STANDARD_NAME]
        E5landata.missing_value=-32767
        E5landata.long_name = __meteo_vars_config[vr][__KEY_LONG_NAME]
        E5landata.units = __meteo_vars_config[vr][__KEY_UNIT]
        #E5landata.valid_min=(int(np.min(rfile))-add_offset)/scale_factor 
        #E5landata.valid_max=(int(np.max(rfile))-add_offset)/scale_factor-100
        
            #E5landata.valid_min=__meteo_vars_config[vr][__KEY_VMIN]
            #E5landata.valid_max=__meteo_vars_config[vr][__KEY_VMAX]   
        
        E5landata.scale_factor=scale_factor        
        E5landata.add_offset=add_offset
        E5landata.set_auto_maskandscale(False)  

        print(E5landata.units)    
        E5landata.grid_mapping='wgs_1984'
        E5landata.esri_pe_string='GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]]'
         

    
        
        #kiko=np.stack(rfile['latitude'])
        #koki=np.stack(rfile['longitude'])
        kaka=len(rfile['time'])
    
        latitude[:]=np.arange(72.25, 22.15, -0.1)    ##-89.999 999 999 999 9716  -89.999 999 999 999 9858 
        longitude[:]=np.arange(-25.25, 50.35, 0.1)
        time[:]=np.arange(kaka)
        
        #wala=np.isnan(results)
        #print(wala)
        #wali=np.asmatrix(results,dtype='uint8')
        __VALUE_NAN=-32767
        resulti=rfile.fillna(__VALUE_NAN)
        #lona=round(__VALUE_NAN  * scale_factor +add_offset,1)
        #resulti[np.isnan(resulti)] = __VALUE_NAN  * scale_factor +add_offset
        #if vr in ['ssr','str']:
            #resulti=resulti.astype('f4')
        #else:
            #resulti=resulti.astype('f4')


   
    
    
    
        for t in range(kaka):
            E5landata[t,:,:]=resulti[t,:,:]
        

    
        ####print(E5landata)
        nf2.close()
        resulti=[]


# Uncomment years as required
years =  [
'1981',
'1982', '1983', '1984',
'1985', '1986', '1987',
'1988', '1989', '1990',
#'1991', '1992', '1993',
#'1994', '1995', '1996',
#'1997', '1998', '1999',
#'2000', '2001', '2002',
'2003', 
'2004', '2005','2006',
'2007', '2008','2009', 
'2010', '2011',
'2012', '2013', '2014',
'2015', '2016', '2017',
'2018', '2019', '2020',
]
months = [ "01",
          "02", "03", "04", "05", "06", 
         "07", "08", "09", "10", "11", "12"
         ] 


# select your variable(s); name must be a valid ERA5 CDS API name 
varconf=['u10', 'v10','t2m','str','ssr','tp','d2m']

# define names of new variables
tvar = ['ws', 'ta','rg','rn','td','tp']
var = tvar

# parameter linked to the way data have been downloaded
# fi=1 // all duration
# fi=2 // 1991 - 2002
# fi=3 // 1981-1990 & 2003-2020
fi=3
for yr in years:

    for mo in months:
            if fi==1:
                hourly_v = xr.open_dataset(os.getcwd() + '/hourly/e5l_'+ yr + '_' + mo + '.nc')
            if fi==2:
                hourly_v = xr.open_dataset(os.getcwd() + '/hourly/e5l_tp-2d_'+ yr + '_' + mo + '.nc') 
            if fi==3:
                hourly_v = xr.open_dataset(os.getcwd() + '/hourly/e5l_tp_'+ yr + '_' + mo + '.nc') 
            
            vax=list(hourly_v.keys())
            print("the variables in this file are " + ' | '.join(vax)) 
            print('month ' + mo )
            start = time.time()
            
            # create mean wind speed from u and v and extract daily mean
            if 'u10' in vax and 'v10' in vax:
                v10= hourly_v['v10']
                u10 = hourly_v['u10']
                wind= wind_uv_to_spd(u10,v10)
                wind=wind.rename({'ws'})
                daily_w = wind.resample(time='D').mean('time')
                if mo=="01":
                    dw=daily_w
                else:
                    dw = concat([dw,daily_w],dim='time')
                    
            # daily accumulation of precipitation   
            if 'tp' in vax:

                tp= hourly_v['tp'] 
                #convertion to mm
                tp=tp*1000
                daily_pr=tp.resample(time='D').max('time')   
                daily_pr=daily_pr.rename({'tp'})
                if mo=="01":
                    dtp=daily_pr
                else:
                    dtp = concat([dtp,daily_pr],dim='time')
                
              # precipitation: calculate sum with frequency of 24h and multiply by 1000
              # precipitation value is for the day before
              
            # daily mean of tempearature
            if 't2m'in vax:
                t2m=hourly_v['t2m']
                #convert Kelvin to degrees C
                t2m=t2m-273.15
                daily_t2m = t2m.resample(time='D').mean('time')
                daily_t2m=daily_t2m.rename({'ta'})
                if mo=="01":
                    dta=daily_t2m
                else:
                    dta = concat([dta,daily_t2m],dim='time')
                    
            # daily means of surface net thermal radiation and surface net solar radiation    
            if 'str'in vax:
                sstr=hourly_v['str']
                daily_rn = sstr.resample(time='D').sum('time',min_count=4) 
                daily_rn=daily_rn.rename({'rn'})
                if mo=="01":
                    drn=daily_rn
                else:
                    drn = concat([drn,daily_rn],dim='time')
                
            if 'ssr'in vax:
                ssr=hourly_v['ssr']
                daily_rg = ssr.resample(time='D').sum('time',min_count=4) 
                daily_rg= daily_rg.rename({'rg'})
                if mo=="01":
                    drg=daily_rg
                else:
                    drg = concat([drg,daily_rg],dim='time')
                
            # daily mean of dew point temperature
            if 'd2m'in vax:
                d2m=hourly_v['d2m']
                #convert Kelvin to degrees C
                d2m=d2m-273.15
                daily_d2m = d2m.resample(time='D').mean('time')
                daily_d2m=daily_d2m.rename({'td'})
                if mo=="01":
                    dtd=daily_d2m
                else:
                    dtd = concat([dtd,daily_d2m],dim='time')
    
            #dailymax_t2m = t2m.resample(time='D').max('time')   
            #dailymax_t2m=dailymax_t2m.rename({'tx'})
            
         
            #dailymin_t2m = t2m.resample(time='D').min('time') 
            #dailymin_t2m=dailymin_t2m.rename({'tn'})

            end=time.time()
            print(end - start)


#%%
    # generate new yearly netcdf of daily values of the variables using scale factor ans add offset
        
    if 'u10' in vax and 'v10' in vax:
        dt=dw
        vr=var[0]
        scale_factor, add_offset = __meteo_vars_config[vr][__KEY_SCALE_FACTOR],__meteo_vars_config[vr][__KEY_OFFSET]
        dtx=np.round((dt - add_offset) / scale_factor)
        print ('Start generating netcdf file for variable: '+ vr)
        print(scale_factor)
        
        writenetcdf(yr,vr,dtx,scale_factor,add_offset)
    if 't2m' in vax:
        dt=dta
        vr=var[1]
        scale_factor, add_offset = __meteo_vars_config[vr][__KEY_SCALE_FACTOR],__meteo_vars_config[vr][__KEY_OFFSET]
        dtx=np.round((dt - add_offset) / scale_factor)
        print ('Start generating netcdf file for variable: '+ vr)
        print(scale_factor)
        
        writenetcdf(yr,vr,dtx,scale_factor,add_offset)
    if 'ssr' in vax:
        dt=drg
        vr=var[2]
        scale_factor, add_offset = __meteo_vars_config[vr][__KEY_SCALE_FACTOR],__meteo_vars_config[vr][__KEY_OFFSET]
        dtx=np.round((dt - add_offset) / scale_factor)
        print ('Start generating netcdf file for variable: '+ vr)
        print(scale_factor)
        
        writenetcdf(yr,vr,dtx,scale_factor,add_offset)
        
    if 'str' in vax:
        dt=drn
        vr=var[3]
        scale_factor, add_offset = __meteo_vars_config[vr][__KEY_SCALE_FACTOR],__meteo_vars_config[vr][__KEY_OFFSET]
        dtx=np.round((dt - add_offset) / scale_factor)
        print ('Start generating netcdf file for variable: '+ vr)
        print(scale_factor)
        
        writenetcdf(yr,vr,dtx,scale_factor,add_offset)
    if 'd2m' in vax:
        dt=dtd
        vr=var[4]
        scale_factor, add_offset = __meteo_vars_config[vr][__KEY_SCALE_FACTOR],__meteo_vars_config[vr][__KEY_OFFSET]
        dtx=np.round((dt - add_offset) / scale_factor)
        print ('Start generating netcdf file for variable: '+ vr)
        print(scale_factor)
        
        writenetcdf(yr,vr,dtx,scale_factor,add_offset)
        
    if 'tp' in vax:
        dt=dtp
        vr=var[5]
        scale_factor, add_offset = __meteo_vars_config[vr][__KEY_SCALE_FACTOR],__meteo_vars_config[vr][__KEY_OFFSET]
        dtx=np.round((dt - add_offset) / scale_factor)
        print ('Start generating netcdf file for variable: '+ vr)
        print(scale_factor)
        
        writenetcdf(yr,vr,dtx,scale_factor,add_offset)