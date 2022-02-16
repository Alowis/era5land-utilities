# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 12:30:25 2021

@author: tilloal
"""
# Import the os module
import os
import rioxarray # for the extension to load
import pandas as pd
import xarray as xr
import numpy as np
import math
from scipy.interpolate import griddata
from xarray import concat
from netCDF4 import Dataset
import netCDF4 as nc
import time

## Script 3/4

# script that allows to reggrid yearly files of daily ERA5* variables to another resolution
# Original grid: 0.1 degrees
# Objective grid: 1 arc minute (approx 0.016667 deg)
 
# Requires:
# 1) yearly files of daily data for the specified variables in the specified years
# (Outputs from "02_ERA5land_yearly_files.py")
 
# Output:
# 1) separated netCDF file for chosen daily variable for each year at the objective resolution


# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))

# key text file containing (1) proxy key, (2) CDS API key, (3) working directory of the data
keys=[]
with open('Keys.txt') as f:
    for line in f:            
        keys.append(line) 
        
# Change the current working directory
os.chdir(keys[2])

# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))

def compute_scale_and_offset(min, max, n):
    # stretch/compress data to the available packed range
    scale_factor = (max - min) / (2 ** n - 1)
    # translate the range to be symmetric about zero
    add_offset = min + 2 ** (n - 1) * scale_factor
    return (scale_factor, add_offset)

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

#%%
# function to write a new netcdf file with objective grid
def regridnetcdf(yr, vr, dt, scale_factor, add_offset,obje, facto):    
    
    tziz=(dt['time'])
    d1=(tziz[0])
    d2=np.datetime_as_string(d1, unit='D')
    tunits= 'days since ' + d2
    fi = nc.Dataset(obje)
    lo=np.asarray(fi.variables['lon'])
    la=np.asarray(fi.variables['lat'])
    d1=os.getcwd() + '\\'
    if facto==True:
        d2=  vr + '\\1_arcmin\e5ld_1arcmin_' + vr + '_'+ yr + '.nc'
    if facto==False:
        d2=  vr + '\\1_arcmin\e5ld_1min_lvap_' + vr + '_'+ yr + '.nc'
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
    aa=dt.shape
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
    E5landata = nf2.createVariable(vr, 'i2',('time','lat','lon',), zlib=True, complevel=4, fill_value=-9999)
    E5landata.standard_name = __meteo_vars_config[vr][__KEY_STANDARD_NAME]
    E5landata.missing_value=-9999
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
    __VALUE_NAN=-9999
    
    #kiko=np.stack(rfile['latitude'])
    #koki=np.stack(rfile['longitude'])
    kaka=len(dt['time']) 
    latitude[:]=la ##-89.999 999 999 999 9716  -89.999 999 999 999 9858 
    #latitude[:]=np.arange(72.25, 22.75, -0.016666667) 
    longitude[:]=lo
    #longitude[:]=np.arange(-25.25, 50.25, 0.016666667)
    
    #daily variable need to be the accumulation of the previous day
    if vr in ['ssr','str', 'ssrd', 'tp']:
        time[:]=np.arange(kaka)
    else:
        time[:]=np.arange(kaka) +1
    #wala=np.isnan(results)
    #print(wala)
    #wali=np.asmatrix(results,dtype='uint8')

    #resulti[np.isnan(resulti)]=(__VALUE_NAN)
    #lona=round(__VALUE_NAN  * scale_factor +add_offset,1)

    #resulti=resulti.astype('i2')
    for t in range(kaka):
        E5landata[t,:,:]=dt[t,:,:]
    
    aa=E5landata
    ####print(E5landata) 
    nf2.close()
    dt=[]

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))


#%%

# compressions factors for variables tp, ws, rg, rn, td
factorz=pd.read_csv('compression_factors.csv')
facto=False


if facto==True:
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
        'rgd' : {__KEY_UNIT : 'J/m2/d', __KEY_STANDARD_NAME : 'ssrd', __KEY_LONG_NAME : 'surface_downward_solar_radiation',
            __KEY_OFFSET :0, __KEY_SCALE_FACTOR : 1000},
        'rn' : {__KEY_UNIT : 'J/m2/d', __KEY_STANDARD_NAME : 'str', __KEY_LONG_NAME : 'surface_net_thermal_radiation',
                __KEY_OFFSET : factorz['add_offset'].values[3], __KEY_SCALE_FACTOR : factorz['scale_factor'].values[3]},
    }

if facto==False:
    __meteo_vars_config = {
        'tn' : {__KEY_UNIT : 'celcius', __KEY_STANDARD_NAME : 'tn', __KEY_LONG_NAME : 'min_temperature',
                __KEY_OFFSET : 0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : -700, __KEY_VMAX : 700},
        'tp' : {__KEY_UNIT : 'mm', __KEY_STANDARD_NAME : 'tp', __KEY_LONG_NAME : 'total_precipitation',
                __KEY_OFFSET : 0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : -1, __KEY_VMAX : 7000},
        'ta' : {__KEY_UNIT : 'celcius', __KEY_STANDARD_NAME : 'ta', __KEY_LONG_NAME : 'mean_temperature',
                __KEY_OFFSET : 0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : -700, __KEY_VMAX : 700},
        'td' : {__KEY_UNIT : 'celcius', __KEY_STANDARD_NAME : 'td', __KEY_LONG_NAME : 'mean_dewpoint_temperature',
                __KEY_OFFSET : 0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : -700, __KEY_VMAX : 700},
        'tx' : {__KEY_UNIT : 'celcius', __KEY_STANDARD_NAME : 'tx', __KEY_LONG_NAME : 'max_temperature',
                __KEY_OFFSET : 0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : -400, __KEY_VMAX : 400},
        'ws' : {__KEY_UNIT : 'm/s', __KEY_STANDARD_NAME : 'ws', __KEY_LONG_NAME : 'avg_wind_speed',
                __KEY_OFFSET : 0, __KEY_SCALE_FACTOR : 0.1, __KEY_VMIN : 0, __KEY_VMAX : 45},
        'rg' : {__KEY_UNIT : 'J/m2/d', __KEY_STANDARD_NAME : 'ssr', __KEY_LONG_NAME : 'surface_net_solar_radiation',
                __KEY_OFFSET : 0, __KEY_SCALE_FACTOR : 1000},
        'rgd' : {__KEY_UNIT : 'J/m2/d', __KEY_STANDARD_NAME : 'ssrd', __KEY_LONG_NAME : 'surface_downward_solar_radiation',
            __KEY_OFFSET :0, __KEY_SCALE_FACTOR : 1000},
        'rn' : {__KEY_UNIT : 'J/m2/d', __KEY_STANDARD_NAME : 'str', __KEY_LONG_NAME : 'surface_net_thermal_radiation',
                __KEY_OFFSET : 0, __KEY_SCALE_FACTOR : 1000},
}
#%%
# Uncomment years as required
years =  [
'1981',
'1982', '1983', '1984',
'1985', '1986', '1987','1988', 
'1989', 
'1990','1991', 
'1992', '1993',
'1994', '1995', '1996',
'1997', 
'1998', '1999','2000',
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

months = [ "01",
          "02", "03", "04", "05", "06", 
         "07", "08", "09", "10", "11", "12"
         ] 

#var = ['ws', 'ta','td','rg','rn','rgd']
#var = ['rgd']
var= ['tp']

# file with the objective grid
obje="Source.nc"

# loop for each year and each variable
for yr in years:
    for vr in var:
        ox=[]
        file_obj = nc.Dataset(vr + "\\0.1_deg\e5ldx_" + vr + "_" + yr+ ".nc")
        #if vr in ['ssr','str']:
        scale_factor, add_offset = __meteo_vars_config[vr][__KEY_SCALE_FACTOR],__meteo_vars_config[vr][__KEY_OFFSET]
        print(file_obj.variables[vr].add_offset)
        print(add_offset)
        
        #scale_factor = file_obj.variables[vr].scale_factor
        Source = xr.open_dataset("Source.nc")
        Target = xr.open_dataset(vr + "\\0.1_deg\e5ldx_" + vr + "_" + yr+ ".nc",mask_and_scale=True)
        #Target = xr.open_dataset(vr + "\\0.1_deg\e5lf_" + vr + "_" + yr+ ".nc")
        obj=Source['area']
        vt=len(Target['time'])
        start = time.time()
        vip=Target[vr]
        vip=np.round((vip - add_offset) / scale_factor)
        vip=vip.rio.write_crs(4326)
        vip = vip.where(vip!=-32767,np.NaN)
        vip=vip.rio.write_nodata('nan')
        
        # nearest neighbors to fill gaps (non continental tiles)
        print("filling gaps")
        ox=vip.rio.interpolate_na(method='nearest').astype('f2')
        #ox=ox.astype('f4')
        #interpolate to new grid with nearest neighbors 
        print('interpolating')
        ox=ox.interp_like(obj, method='nearest')
        ox=ox.astype('i2')
        #ox.plot()
        # condition to draw the new coastline according to the objective grid
        condition = obj.notnull()
        # fill the values with NA where condition is false
        ox = ox.where(condition, -9999)
        end=time.time()
        print(end - start)

        print ('Start generating netcdf file for variable: '+ vr)
        regridnetcdf(yr,vr,ox,scale_factor,add_offset,obje, facto)

