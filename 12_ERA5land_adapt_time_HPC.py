# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 11:30:03 2021

@author: tilloal
"""

from netCDF4 import Dataset
# Import the os module
import os



## 12. Extra script (ONLY WORKS ON LINUX)

# script that shifts the variable dates by one 
# Useful for inputs of LISVAP and LISFLOOD
 
# Requires:
# 1)  files of of daily data for the specified variables in the specified years
 
# Output:
# 1) netCDF file with shifted date by one day


# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))

# Change the current working directory
os.chdir('/home/tilloal/lisflood-lisvap/')
#%%
vr='ta'
yr='2003'
d1=os.getcwd() + '/'
d2=  'e5ld_1arcmin_' + vr + '_'+ yr + '.nc'
#d3=  vr + '\\1_arcmin\e5lv_1arcmin_' + vr + '_'+ yr + '.nc'
namenc=d1+d2
#nemo=d1+d3


nf2=Dataset(namenc,mode="a",format='NETCDF4_CLASSIC')



tt1=nf2.variables['time'][:]



tt1=tt1+1
nf2.variables['time']=tt1
print(nf2.variables['time'][:])
#%%
nf2.close()