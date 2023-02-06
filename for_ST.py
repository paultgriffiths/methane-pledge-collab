#!/usr/bin/env python
# coding: utf-8

import regionmask
import xarray as xr
import geopandas

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pyproj

###
# sanity checks
###

try: 
    geopandas.__version__ == 0.8
except:
    print('for some reason geopandas 0.8 works on JASMIN, later versions segfault')
    

try: 
    regionmask.__version__ == 0.9
except:
    print('we need regionmask 0.9 ')

### function defs
def area_grid(lat,lon):
    '''
    Function to calculate surface area per gridbox
    Units: m2
    S = R^2*(lon2-lon1)*(sin lat2 - sin lat1)
    lon in radians, R = 6371 km
    '''
    import numpy as np
    Pi           = np.float128(3.141592653589793238462643383279)
    Earth_Radius = np.float128(6371.0*1.0E3)#equator radius:6378.1*1E3
    lat_bound    = np.float128(89.999999999999999999999999)
    lon          = np.float128(lon)
    lat          = np.float128(lat)
    rlon         = (lon[:]/np.float128(180.0))*Pi
    rlat         = (lat[:]/np.float128(180.0))*Pi
    dlat         = (rlat[1] - rlat[0])/2.0
    dlon         = (rlon[1] - rlon[0])/2.0
    #
    area = np.zeros((len(rlat),len(rlon)),np.float128)
    j=0
    while j < len(rlat):
        if (lat[j] >= lat_bound):
            lat1 = rlat[j]
            lat2 = rlat[j] - dlat/2.0
        elif (lat[j] <= -1.0*lat_bound):
            lat1 = rlat[j] + dlat/2.0
            lat2 = rlat[j]
        else:
            lat1 = rlat[j] + dlat
            lat2 = rlat[j] - dlat
        i=0
        while i < len(rlon):
            lon1 = rlon[i] - dlon
            lon2 = rlon[i] + dlon
            area[j,i] = (Earth_Radius**2)*(abs(np.sin(lat1)-np.sin(lat2))*abs(lon1-lon2))
            i += 1
        j += 1
    return area

# JASPY bespoke kernel doesn't pick up on the correct pyproj dir...
try:
    pyproj.datadir.set_data_dir('/home/users/ptg21/anaconda3/envs/geopandas-0.8_2022-10-18/share/proj')
except:
    print('please set pyproj directory as above, will need to change USER and ENVIRONMENT name')

# find lat/lons of China etc 
CN_index =  139# regionmask.defined_regions.natural_earth_v4_1_0.countries_110.map_keys("CN")
RU_index =  regionmask.defined_regions.natural_earth_v4_1_0.countries_110.map_keys("RUS")
IN_index =  regionmask.defined_regions.natural_earth_v4_1_0.countries_110.map_keys("IND")
PL_index =  regionmask.defined_regions.natural_earth_v4_1_0.countries_110.map_keys("PL")
AUS_index = regionmask.defined_regions.natural_earth_v4_1_0.countries_110.map_keys("AU")
TKM_index = regionmask.defined_regions.natural_earth_v4_1_0.countries_110.map_keys("TM")
IRN_index = regionmask.defined_regions.natural_earth_v4_1_0.countries_110.map_keys("IRN")


##
## create an array to work on/with that will hold the emissions in kgm-2s-1
## this will be the xarray dataset that gets written to disk with the new emms
##
CMIP6_ems =  xr.load_dataset("/gws/nopw/j04/htap2/ptg21/methane-pledge-emissions/CH4_anthropogenic_2014_2101_time_series.nc")
MP_ems = CMIP6_ems.copy(deep=True)

##
## create an xarray dataset that will store the mask
##
mask_for_ems = xr.load_dataset("/gws/nopw/j04/htap2/ptg21/methane-pledge-emissions/CH4_anthropogenic_2014_2101_time_series.nc")
mask_for_ems['emissions_mask_CH4']= mask_for_ems['emissions_CH4'].copy(deep=True)

##
## make a mask - modify this to suit
## at present the mask is time-varying
## 
mask_for_ems.emissions_mask_CH4.data[...]=0.

##
## create time-varying mask
## 1. for first 9 years do nothing
##
mask_for_ems['emissions_mask_CH4'][0:108,...]=1.

##
## 2. after this decrease linearly to 0.7 everywhere over 7 years
##
for ivar in range(0,8):
    mask_for_ems['emissions_mask_CH4'][108+(ivar*12):108+(ivar+1)*12,...]=1.-(0.3/7.)*ivar

##
## 3. after this stay at 0.7
##
mask_for_ems['emissions_mask_CH4'][204:,...]=0.7

## now apply country masking to emissions
## country mask from UKCA grid
countrymask = regionmask.defined_regions.natural_earth_v4_1_0.countries_110.mask(mask_for_ems, lon_name='longitude', lat_name='latitude')

##
## the holdouts do nothing so they have a mask = 1.0 for all time
## overwrite existing mask values with 1.0 everywhere
##
mask_for_ems = mask_for_ems.where(
    (
    (countrymask != 139)&
    (countrymask != RU_index) &
    (countrymask != IN_index) & 
    (countrymask != PL_index) & 
    (countrymask != AUS_index) &
    (countrymask != TKM_index) &
    (countrymask != IRN_index) 
    )
    ,  other=1.)

##
## extract mask values from the xarray dataset- needs to be as an np.array as will use it below
##
mask_for_ems_np = np.array(mask_for_ems.emissions_mask_CH4.data)

# construct emissions as emissions * mask
MP_ems['emissions_CH4'].data[:] = mask_for_ems_np*CMIP6_ems.emissions_CH4.data
MP_ems['emissions_CH4'].compute()

# write to netcdf
MP_ems.to_netcdf('CH4_anthropogenic_2014_2101_time_series_methane_pledge_CN_RU_IN_AUS_PL_IR_TK.nc')

## test plots - set make_test_plots True for sanity check
make_test_plots = True

if make_test_plots: 
    plt.figure(figsize=(12,4),dpi=200)
    plt.subplot(1,3,1)
    plt.pcolormesh(mask_for_ems['emissions_mask_CH4'][0,0,...], vmin=0.0, vmax=1.)
    plt.colorbar()
    plt.title('mask at timestep 0')

    plt.subplot(1,3,2)
    plt.pcolormesh(mask_for_ems['emissions_mask_CH4'][225,0,...], vmin=0.0, vmax=1.)
    plt.colorbar()
    plt.title('mask at timestep 225')
    
    plt.subplot(1,3,3)
    plt.pcolormesh(mask_for_ems['emissions_mask_CH4'][-1,0,...], vmin=0.0, vmax=1.) 

    plt.colorbar()
    plt.title('mask at final timestep')
    plt.savefig('test1.png')

if make_test_plots:
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    f, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))#, figsize=(20,20))
    ax.coastlines()

    regionmask.defined_regions.natural_earth_v4_1_0.countries_110.plot(
        ax=ax, add_label=False, line_kws=dict(lw=0.5, color="0.5"))
    mask_for_ems.isel(time=150).emissions_mask_CH4[-1,:,:].plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(),vmin=0.7, vmax=1.)
    plt.title('ignore units - is just a mask')
    plt.savefig('test2.png')

