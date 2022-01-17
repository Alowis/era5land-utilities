# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 10:06:33 2021

@author: tilloal
"""

import netCDF4 as nc
import xarray as xr
from netCDF4 import Dataset
import os
import numpy as np
import pandas as pd

#%%
## Script 4/4

# script that concatenates yearly files of daily ERA5* variables into 
# a single file for the whole specified
 
# Requires:
# 1) monthly files of daily data for the specified variables in the specified years
# (Outputs from "03_ERA5land_interpolate.py")
 
# Output:
# 1) separate netCDF file for chosen daily variable for the specified period


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


# key text file containing (1) proxy key, (2) CDS API key, (3) working directory of the data
keys=[]
with open('Keys.txt') as f:
    for line in f:            
        keys.append(line) 
        
# Change the current working directory
os.chdir(keys[2])

# Print the current working directory

print("Current working directory: {0}".format(os.getcwd()))

# compressions factors for variables tp, ws, rg, rn, td
factorz=pd.read_csv('compression_factors.csv')

# file with objective long lat grid
src="Source.nc"

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

var = ['ws', 'ta','td','rg','rn','tp']
# select variable for which a file will be created (from var above)
vr="td"

# Uncomment years as required
years =  [
#'1981',
'1982', '1983', '1984',
'1985', '1986', '1987',
'1988', '1989', '1990',
'1991', 
'1992', '1993',
'1994', '1995', '1996',
'1997', '1998', '1999','2000',
'2001', '2002',
'2003', '2004', '2005',
'2006', '2007', '2008',
'2009',
'2010', '2011',
'2012',
'2013', '2014',
'2015', '2016', '2017',
'2018', '2019', '2020',

]

scale_factor, add_offset = __meteo_vars_config[vr][__KEY_SCALE_FACTOR],__meteo_vars_config[vr][__KEY_OFFSET]

ft=[]
for yr in years:
    d1=os.getcwd() + '/' + vr + '/1_arcmin/'
    d1l=os.getcwd() + '/' + vr + '/0.1_deg/'
    d2=  'e5ld_1arcmin_' + vr + '_'+ yr + '.nc'
    d2l=  'e5ldx_' + vr + '_'+ yr + '.nc'
    f=d1+d2
    fl=d1l+d2l
    ft.append(f)

ds = xr.open_mfdataset(ft,combine = 'nested', concat_dim="time",mask_and_scale=False)
dsx=ds[vr]

#%%
# tweaked writenetcdf function for large files
def writeMega(vr, dsx, scale_factor, add_offset,src):   
        #rfile = dsx.copy()
        tziz=(dsx['time'])
        d1=(tziz[0])

        fi = nc.Dataset(src)
        lo=np.asarray(fi.variables['lon'])
        la=np.asarray(fi.variables['lat'])

        d2=np.datetime_as_string(d1, unit='D')
        d3=" 00:00:00"
        d4=d2+d3
        tunits= 'days since ' + d2
        d1=os.getcwd() + '\\'
        d2=  vr + '\\1_arcmin\e5l_1amin_1981_2020_' + vr + '.nc'
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
        aa=dsx.shape
        nf2.createDimension('lon', aa[2])
        nf2.createDimension('lat', aa[1])
        nf2.createDimension('time', None)
         #Variables
        longitude = nf2.createVariable('lon','f8',('lon',), complevel=4, zlib=True)  ###('lon','f8',('lon',))
        longitude.standard_name= 'Longitude'
        longitude.long_name= 'Longitude'
        longitude.units =  'degrees_east'
        
        latitude = nf2.createVariable('lat','f8',('lat',), complevel=4, zlib=True)  ###('lon','f8',('lon',))
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
        E5landata = nf2.createVariable(vr, 'i2',('time','lat','lon',), zlib=True, complevel=4, fill_value=-32767)
        E5landata.standard_name = __meteo_vars_config[vr][__KEY_STANDARD_NAME]
        E5landata.missing_value=-32767
        E5landata.long_name = __meteo_vars_config[vr][__KEY_LONG_NAME]
        E5landata.units = __meteo_vars_config[vr][__KEY_UNIT]
        #E5landata.valid_min=(int(np.min(rfile))-add_offset)/scale_factor 
        #E5landata.valid_max=(int(np.max(rfile))-add_offset)/scale_factor-100       
        E5landata.scale_factor=scale_factor        
        E5landata.add_offset=add_offset
        print(E5landata.units)    
        E5landata.set_auto_maskandscale(False)  
        E5landata.grid_mapping='wgs_1984'
        E5landata.esri_pe_string='GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]]' 
        __VALUE_NAN=-32767
        
        #kiko=np.stack(rfile['latitude'])
        #koki=np.stack(rfile['longitude'])
        kaka=len(dsx['time']) 
        latitude[:]=la ##-89.999 999 999 999 9716  -89.999 999 999 999 9858 
        #latitude[:]=np.arange(72.25, 22.75, -0.016666667) 
        longitude[:]=lo
        #longitude[:]=np.arange(-25.25, 50.25, 0.016666667)
        time[:]=np.arange(kaka)
        #wala=np.isnan(results)
        #print(wala)
        #wali=np.asmatrix(results,dtype='uint8')
    
        #resulti[np.isnan(resulti)]=(__VALUE_NAN)
        #lona=round(__VALUE_NAN  * scale_factor +add_offset,1)

        #resulti=resulti.astype('i2')
        for t in range(kaka):
            E5landata[t,:,:]=dsx[t,:,:]
     
        aa=E5landata
        ####print(E5landata) 
        nf2.close()
        dsx=[]
#%%
print ('Start generating netcdf file for variable: '+ vr)
writeMega(vr,dsx,scale_factor,add_offset,src)