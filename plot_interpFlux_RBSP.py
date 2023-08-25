#%%
import json
import pytplot
import pytz
import pandas as pd

from pysatdata.loaders.load import *
from pysatdata.utils.interpolate_flux_rbsp import *
from pysatdata.utils.plotFunc.plot_funtions import *


# File with the data sources and local destinations
config_file_sat = 'config_file.json'

trange_plot=['2016-07-06', '2016-07-10'] # time range for ploting
trange=['2016-07-05', '2016-07-11'] # time range for interpolated data (requires a little bit more to overcome the edge problemns in interpolation)
#%%
# test Download MAgEmphemeri data
paramLoadSat = {"satellite": 'rbsp', "probe": 'a', "level": 'TS04', "rel": "rel03",
                "instrument": 'MagEphem', "datatype": 'definitive'}

varss_ephem = load_sat(trange=trange, satellite=paramLoadSat["satellite"],
                     probe=[paramLoadSat["probe"]], level=paramLoadSat["level"], rel=paramLoadSat["rel"],
                     instrument=paramLoadSat["instrument"],datatype=paramLoadSat["datatype"], downloadonly=True,
                     testRemotePath=True,
                     usePandas=False, usePyTplot=True)
#%%
# test loading merged Ace data
# # define the details for OMNI datsa (Only if you need it)
# paramLoadSat = {"satellite": 'ace', "probe": 'ace',"instrument": 'merged', "datatype": '4_min'}
# pytplot.del_data()
# varss_aceSwe = load_sat(trange=trange, satellite=paramLoadSat['satellite'],
#                          probe=[paramLoadSat['probe']], rel='rel03',
#                          instrument=paramLoadSat['instrument'],datatype=paramLoadSat['datatype'],
#                          config_file=config_file_sat, downloadonly=False,
#                          usePandas=True, usePyTplot=False)

#%%

paramLoadSat = {"satellite": 'rbsp', "probe": 'a', "level": '3', "rel": "rel03",
                "instrument": 'rept', "datatype": 'sectors'}

# The command below will download the data in the "sat_data" repository in the root directory of
# your computer (This is the default, you can change it in the config_file.json file). 
# The variables will be stored as pytplot format, the same way as pyspedas does.
varss_rept = load_sat(trange=trange, satellite=paramLoadSat["satellite"],
                     probe=[paramLoadSat["probe"]], level=paramLoadSat["level"], rel=paramLoadSat["rel"],
                     instrument=paramLoadSat["instrument"],datatype=paramLoadSat["datatype"],
                     config_file=config_file_sat, downloadonly=False,
                     usePandas=False, usePyTplot=True)
#%%
# Since the variables are stored in memory we need to open to perform
# the cleaning and select the desired chanel
quants_fedu_rept = pytplot.data_quants['FEDU'] # FEDU is the high energy flux electron of Rept instrument
flux_rept = quants_fedu_rept.values # flux_rept is a 3D variavle (time x pitchAngle x energy)
flux_rept[flux_rept == -9999999848243207295109594873856.000] = np.nan
flux_rept[flux_rept == -1e31] = np.nan

# Get the mean value over all pitch angles
flux_rept_spec = np.nanmean(flux_rept, axis=1)
l_rept = pytplot.data_quants['L_star'].values 
l_rept[l_rept == -9999999848243207295109594873856.000] = np.nan
l_rept[l_rept == -1e31] = np.nan
time_dt_rept = quants_fedu_rept.coords['time'].values
# time_dt_rept = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time_rept]
spec = quants_fedu_rept.coords['spec_bins'].values
v1 = quants_fedu_rept.coords['v1'].values

#%%
# define the details for OMNI datsa (Only if you need it)
paramLoadSat = {"satellite": 'omni', "probe": 'omni',"instrument": 'omni_cdaweb', "datatype": 'hro_1min'}
pytplot.del_data()
varss_aceSwe = load_sat(trange=trange, satellite=paramLoadSat['satellite'],
                         probe=[paramLoadSat['probe']], rel='rel03',
                         instrument=paramLoadSat['instrument'],datatype=paramLoadSat['datatype'],
                         config_file=config_file_sat, downloadonly=False,
                         usePandas=False, usePyTplot=True)

#%%
quants_Swe = pytplot.data_quants['proton_density']
time_dt_swe = quants_Swe.coords['time'].values
# time_dt_swe = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time_swe]
nP = pytplot.data_quants['proton_density'].values
bgse_x = pytplot.data_quants['BX_GSE'].values
bgse_y = pytplot.data_quants['BY_GSM'].values
bgse_z = pytplot.data_quants['BZ_GSM'].values
b_total = np.sqrt((bgse_x**2 + bgse_y**2 + bgse_z**2))
flow_speed = pytplot.data_quants['flow_speed'].values
imf = pytplot.data_quants['IMF'].values
ae_index = pytplot.data_quants['AE_INDEX'].values
symH_index = pytplot.data_quants['SYM_H'].values

#%%
fluxEnergyChanel = 1 # Deine the energy level to do the interpolation
# Function to interpolate the electron flux and L-shell
x_interp, y_interp, interpFlux = interpolateFluxRbsp(flux_rept_spec[:,fluxEnergyChanel], l_rept, pd.to_datetime(time_dt_rept))
cutLshell = 5.
energyRange = range(0,4)

# Those are the input to plot a graph with the interpolated flux and 
# the solar wind paramenters. You can find the plot funtion at
# pysatdata/utils/plotFunc/plot_funtions.py 
# this is just an example. The most important is to have the 
# interFlux variable and the x, and y with the same resolution

plotFluxParamsDict = {
                      'specEnergy': spec,
                      'time_dt_rept': time_dt_rept,
                      'l_probe': l_rept,
                      'cutLshell': cutLshell,
                      'energyRange': energyRange,
                      'trangeXlim': trange_plot,
                      'fluxEnergyChanel': fluxEnergyChanel,
                      'flux_rept_spec': flux_rept_spec,
                      'xax': x_interp,
                      'yax': y_interp,
                      'maskflux': interpFlux,
                      'time_dt_swe': time_dt_swe,
                      'nP': nP,
                      'bgse_x': bgse_x,
                      'bgse_y': bgse_y,
                      'bgse_z': bgse_z,
                      'b_total': b_total,
                      'flow_speed': flow_speed,
                      'imf': imf,
                      'ae_index': ae_index,
                      'symH_index': symH_index,
                      'fontsize': 16,
                      'probe': 'A',
                      'level': 2,
                      'plotDir': 'rbsp/rept/',
                      'figureIdentify': None,
                      'interpolatedFlux': True,
                      'vmin': 1,
                      'vmax': 1e6,
                      'lshellOrlStar': 'L Shell'
                     }


#%%

plot_classicFluxSWparams(**plotFluxParamsDict)

# %%
