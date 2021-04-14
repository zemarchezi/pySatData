import pandas as pd
from pysatdata.utils.netcdf2tplot import *
from netCDF4 import Dataset
import xarray as xr

#%%

def readData_goes(files, usePyTplot, usePandas, suffix, time):

    for file in files:
        tvars = netcdf_to_tplot(file, suffix=suffix, merge=True, time=time)


    return tvars
