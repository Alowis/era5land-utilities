# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 19:25:05 2021

@author: tilloal
"""
## Extra script

# Uses calculated add offset and scale factor from "02_ERA5land_yearly_files.py"
# to estimate the optimal parameters for the whole period (1981-2020)
 
 
# Output:
# 1) Valid compression parameters for each variable on thw whole period

import netCDF4 as nc
import xarray as xr
from netCDF4 import Dataset
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
# Change the current working directory
os.chdir('Z:/ClimateRun3/ERA5-land//')

# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))

years =  [
'1981',
'1982', '1983', '1984',
'1985', '1986', '1987',
'1988', '1989', '1990',
'1991', '1992', '1993',
'1994', '1995', '1996',
'1997', '1998', '1999',
'2000', '2001', '2002',
#'2003', 
#'2004', '2005','2006',
#'2007', '2008','2009', 
#'2010', '2011',
#'2012', '2013', '2014',
#'2015', '2016', '2017',
#'2018', '2019', '2020'
]


# Initialisation of the relust table
fdt = pd.DataFrame(np.full([1,3], None),
                   columns=['var', 'add_offset', 'scale_factor'])
var=['ws', 'ta','rg','rn','td']
vrd=['tp']

fdt['var']=vrd
i=-1
f, axarr = plt.subplots(5, sharex=True)
f.suptitle('scale factor')

for vr in vrd:
    i=i+1
    adof=[]
    scal=[]
    mev=[]
    mav=[]
    for yr in years:
            file_obj = nc.Dataset(vr + "\\0.1_deg\e5ldx_" + vr + "_" + yr+ ".nc")
            yearly_v = xr.open_dataset(vr + "\\0.1_deg\e5ldx_" + vr + "_" + yr+ ".nc")[vr]
            
            add_offset = file_obj.variables[vr].add_offset
            scale_factor = file_obj.variables[vr].scale_factor
            am=np.mean(yearly_v)
            #ma=np.max(yearly_v)
            
            mev += [am]
            #mav += [ma]
            adof += [add_offset]
            scal += [scale_factor]
            
    #plt.scatter(range(40),scal)
    mep=np.array(mev)
    axarr[i].plot(range(22), mep)
    if vr in ['td']:
        bip=np.argmax(scal[0:25])
        
    if vr in ['tp']:
        bip=np.argmax(scal[0:21])
    else:
        bip=np.argmax(scal)

    
    adf=adof[bip]
    scf=scal[bip]
    fdt['add_offset'][i]=adf
    fdt['scale_factor'][i]=scf
#%%   
j=-1
for ax in axarr.flat:
    j=j+1
    ax.set(ylabel=var[j])
#%%
# export of the reslut in csv format
fdt.to_csv('compression_factors.csv',index=False)