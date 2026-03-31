#%%
import json
import pytplot
import pytz
import pandas as pd

from pysatdata.loaders.load import *
from pysatdata.utils.interpolate_flux_rbsp import *
from pysatdata.utils.plotFunc.plot_funtions import *


trange_plot=['2015-01-02', '2015-01-03'] # time range for ploting
trange=['2015-01-01', '2015-01-04'] # time range for interpolated data (requires a little bit more to overcome the edge problemns in interpolation)
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


paramLoadSat = {"satellite": 'rbsp', "probe": 'b', "level": '3', "rel": "rel03",
                "instrument": 'rept', "datatype": 'sectors'}

# The command below will download the data in the "sat_data" repository in the root directory of
# your computer (This is the default, you can change it in the config_file.json file). 
# The variables will be stored as pytplot format, the same way as pyspedas does.
varss_rept = load_sat(trange=trange, satellite=paramLoadSat["satellite"],
                     probe=[paramLoadSat["probe"]], level=paramLoadSat["level"], rel=paramLoadSat["rel"],
                     instrument=paramLoadSat["instrument"],datatype=paramLoadSat["datatype"],
                     downloadonly=False,
                     usePandas=False, usePyTplot=True)
#%%
# Since the variables are stored in memory we need to open to perform
# the cleaning and select the desired chanel
quants_fedu_rept = pyspedas.get('FEDU') # FEDU is the high energy flux electron of Rept instrument
flux_rept = np.asarray(quants_fedu_rept.y) # flux_rept is a 3D variavle (time x pitchAngle x energy)
flux_rept[flux_rept == -9999999848243207295109594873856.000] = np.nan
flux_rept[flux_rept == -1e31] = np.nan

# Get the mean value over all pitch angles
flux_rept_spec = np.nanmean(flux_rept, axis=1)
l_rept = np.asarray(pyspedas.get('L_star').y)
l_rept[l_rept == -9999999848243207295109594873856.000] = np.nan
l_rept[l_rept == -1e31] = np.nan
time_dt_rept = np.asarray(quants_fedu_rept.times)
# time_dt_rept = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time_rept]
spec = np.asarray(quants_fedu_rept.v2)
v1 = np.asarray(quants_fedu_rept.v1)

#%%
# define the details for OMNI datsa (Only if you need it)
paramLoadSat = {"satellite": 'omni', "probe": 'omni',"instrument": 'omni_cdaweb', "datatype": 'hro_1min'}
pytplot.del_data()
varss_aceSwe = load_sat(trange=trange, satellite=paramLoadSat['satellite'],
                         probe=[paramLoadSat['probe']], rel='rel03',
                         instrument=paramLoadSat['instrument'],datatype=paramLoadSat['datatype'],
                         downloadonly=False,
                         usePandas=False, usePyTplot=True)

#%%
quants_Swe = pyspedas.get('proton_density')
time_dt_swe = np.asarray(quants_Swe.times)
# time_dt_swe = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time_swe]
nP = np.asarray(pyspedas.get('proton_density').y)
bgse_x = np.asarray(pyspedas.get('BX_GSE').y)
bgse_y = np.asarray(pyspedas.get('BY_GSM').y)
bgse_z = np.asarray(pyspedas.get('BZ_GSM').y)
b_total = np.sqrt((bgse_x**2 + bgse_y**2 + bgse_z**2))
flow_speed = np.asarray(pyspedas.get('flow_speed').y)
imf = np.asarray(pyspedas.get('IMF').y)
ae_index = np.asarray(pyspedas.get('AE_INDEX').y)
symH_index = np.asarray(pyspedas.get('SYM_H').y)

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
                      'specEnergy': spec, #eixo y plot do fluxo interpolado
                      'l_probe': l_rept, # l shell ou l* 
                      'cutLshell': cutLshell, # corte do fluxo
                      'energyRange': energyRange, # faixa de energia para plotar os cortes
                      'trangeXlim': trange_plot, # limite do eixo X para o plot [min, max]
                      'fluxEnergyChanel': fluxEnergyChanel, # 
                      'flux_rept_spec': flux_rept_spec, # fluxo de eletrons do rept
                      'time_dt_rept': time_dt_rept, # eixo x para os dados do rept
                      'xax': x_interp, # eixo x para o fluxo interpolado
                      'yax': y_interp, # eixo y para o fluxo interpolado
                      'maskflux': interpFlux, # fluxo interpolado
                      'time_dt_swe': time_dt_swe, # eixo x para os dados do ace
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