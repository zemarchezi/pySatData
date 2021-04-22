import json
import pytplot
import pytz

from pysatdata.loaders.load import *
from pysatdata.utils.interpolate_flux_rbsp import *
from pysatdata.utils.plotFunc.plot_funtions import *

with open('./pysatdata/utils/resources/config_file.json', 'r') as f:
    config_file_sat = json.load(f)

trange0=['2016-03-06', '2016-03-10']
trange=['2016-03-05', '2016-03-11']

#%%
varss_rept = load_sat(trange=trange, satellite='rbsp',
                     probe=['a'], level='l2', rel='rel03',
                     instrument='rept',datatype='h2',
                     config_file=config_file_sat, downloadonly=False,
                     usePandas=False, usePyTplot=True)
quants_fedu_rept = pytplot.data_quants['FEDU']
flux_rept = quants_fedu_rept.values
flux_rept[flux_rept == -9999999848243207295109594873856.000] = np.nan
flux_rept[flux_rept == -1e31] = np.nan

flux_rept_spec = np.nanmean(flux_rept, axis=1)
l_rept = pytplot.data_quants['L_star'].values
l_rept[l_rept == -9999999848243207295109594873856.000] = np.nan
l_rept[l_rept == -1e31] = np.nan
time_rept = quants_fedu_rept.coords['time'].values
time_dt_rept = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time_rept]
spec = quants_fedu_rept.coords['spec_bins'].values
v1 = quants_fedu_rept.coords['v1'].values


pytplot.del_data()
varss_aceSwe = load_sat(trange=trange, satellite='omni',
                     probe=[''], level='hro', rel='rel03',
                     instrument='swe',datatype='1min',
                     config_file=config_file_sat, downloadonly=False,
                     usePandas=False, usePyTplot=True)

quants_Swe = pytplot.data_quants['proton_density']
time_swe = quants_Swe.coords['time'].values
time_dt_swe = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time_swe]
nP = pytplot.data_quants['proton_density'].values
bgse_x = pytplot.data_quants['BX_GSE'].values
bgse_y = pytplot.data_quants['BY_GSE'].values
bgse_z = pytplot.data_quants['BZ_GSE'].values
b_total = np.sqrt((bgse_x**2 + bgse_y**2 + bgse_z**2))
flow_speed = pytplot.data_quants['flow_speed'].values
imf = pytplot.data_quants['IMF'].values
ae_index = pytplot.data_quants['AE_INDEX'].values
symH_index = pytplot.data_quants['SYM_H'].values

#%%
fluxEnergyChanel = 1
xax, yax, maskflux = interpolateFluxRbsp(flux_rept_spec[:,fluxEnergyChanel], l_rept, time_dt_rept)
cutLshell = 5.
energyRange = range(0,4)
plotFluxParamsDict = {
                      'specEnergy': spec,
                      'time_dt_rept': time_dt_rept,
                      'l_probe': l_rept,
                      'cutLshell': cutLshell,
                      'energyRange': energyRange,
                      'trangeXlim': trange0,
                      'fluxEnergyChanel': fluxEnergyChanel,
                      'flux_rept_spec': flux_rept_spec,
                      'xax': xax,
                      'yax': yax,
                      'maskflux': maskflux,
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
                      'plotDir': 'rbsp/rept/'
                     }


#%%

plot_classicFluxSWparams(**plotFluxParamsDict)

