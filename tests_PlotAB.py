#%%
import numpy as np
import pandas as pd
import pytz
from pysatdata.loaders.load import *
from pysatdata.utils.interpolate_flux_rbsp import *
from pysatdata.utils.plotFunc.plot_funtions import *
import datetime
import gc
from pysatdata.utils.library_functions import *
#%%
# dataT = pd.read_csv('/home/jose/python_projects/rbsp_flux_sw_sea/dataEpoch/NewallEpochs_HSS_corrected_2018.csv',
                    # index_col=0)
#
# for n, dd in enumerate(dataT.index):

stringInstant = '2014-04-07'

# stringInstant = dataT.index[0].split(" ")[0]
instDate = datetime.datetime.strptime(stringInstant, '%Y-%m-%d')


inidate = instDate - datetime.timedelta(days = 6)
enddate = instDate + datetime.timedelta(days = 6)
trange0 = [inidate.strftime('%Y-%m-%d'), enddate.strftime('%Y-%m-%d')]
trange= [(inidate - datetime.timedelta(days = 1)).strftime('%Y-%m-%d'),
         (enddate + datetime.timedelta(days = 1)).strftime('%Y-%m-%d')]

config_file_sat = '/home/jose/python_projects/pySatData/pysatdata/resources/config_file.json'



lshellOrlStar = 'L-Star'

paramLoadSat_a = {"satellite": 'rbsp', "probe": 'a', "level": '3', "rel": "rel03",
                "instrument": 'rept', "datatype": 'sectors'}
#%%
varss_rept_a = load_sat(trange=trange, satellite=paramLoadSat_a['satellite'],
                         probe=[paramLoadSat_a['probe']], level=paramLoadSat_a['level'], rel='rel03',
                         instrument=paramLoadSat_a['instrument'],datatype=paramLoadSat_a['datatype'],
                         config_file=config_file_sat, downloadonly=False,
                         usePandas=False, usePyTplot=True)
# varss_rept = load_sat(trange=trange, satellite='rbsp',
#                      probe=[probe], level=level, rel='rel03',
#                      instrument='rept',datatype='sectors',
#                      config_file='./pysatdata/resources/config_file.json', downloadonly=False,
#                      usePandas=False, usePyTplot=True)
quants_fedu_rept_a = pytplot.data_quants['FEDU']
flux_rept_a = quants_fedu_rept_a.values
flux_rept_a[flux_rept_a == -9999999848243207295109594873856.000] = np.nan
flux_rept_a[flux_rept_a == -1e31] = np.nan

flux_rept_spec_a = np.nanmean(flux_rept_a, axis=1)
if lshellOrlStar == 'L-Shell':
    l_rept_a = pytplot.data_quants['L'].values
else:
    l_rept_a = pytplot.data_quants['L_star'].values
l_rept_a[l_rept_a == -9999999848243207295109594873856.000] = np.nan
l_rept_a[l_rept_a == -1e31] = np.nan
time_rept_a = quants_fedu_rept_a.coords['time'].values
time_dt_rept_a = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time_rept_a]
spec_a = quants_fedu_rept_a.coords['spec_bins'].values
v1_a = quants_fedu_rept_a.coords['v1'].values

pytplot.del_data()
paramLoadSat_b = {"satellite": 'rbsp', "probe": 'b', "level": '3', "rel": "rel03",
                "instrument": 'rept', "datatype": 'sectors'}
varss_rept_b = load_sat(trange=trange, satellite=paramLoadSat_b['satellite'],
                         probe=[paramLoadSat_b['probe']], level=paramLoadSat_b['level'], rel='rel03',
                         instrument=paramLoadSat_b['instrument'],datatype=paramLoadSat_b['datatype'],
                         config_file=config_file_sat, downloadonly=False,
                         usePandas=False, usePyTplot=True)
quants_fedu_rept_b = pytplot.data_quants['FEDU']
flux_rept_b = quants_fedu_rept_b.values
flux_rept_b[flux_rept_b == -9999999848243207295109594873856.000] = np.nan
flux_rept_b[flux_rept_b == -1e31] = np.nan

flux_rept_spec_b = np.nanmean(flux_rept_b, axis=1)
if lshellOrlStar == 'L-Shell':
    l_rept_b = pytplot.data_quants['L'].values
else:
    l_rept_b = pytplot.data_quants['L_star'].values
l_rept_b[l_rept_b == -9999999848243207295109594873856.000] = np.nan
l_rept_b[l_rept_b == -1e31] = np.nan
time_rept_b = quants_fedu_rept_b.coords['time'].values
time_dt_rept_b = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time_rept_b]
spec_b = quants_fedu_rept_b.coords['spec_bins'].values
v1_b = quants_fedu_rept_b.coords['v1'].values
#%%
fluxEnergyChanel = 1
xax, yax, maskflux = interpolateFluxRbsp(flux_rept_spec_a[:,fluxEnergyChanel], l_rept_a, time_dt_rept_a)
cutLshell = 5.
energyRange = range(0,4)

energChanel = f'{spec_a[fluxEnergyChanel]}'

plotFluxParamsDict = {'specEnergy': spec_a,
                      'time_dt_rept': [time_dt_rept_a, time_dt_rept_b],
                      'l_probe': [l_rept_a, l_rept_b],
                      'cutLshell': cutLshell,
                      'energyRange': energyRange,
                      'trangeXlim': trange0,
                      'fluxEnergyChanel': fluxEnergyChanel,
                      'fluxAB': [flux_rept_spec_a[:,fluxEnergyChanel], flux_rept_spec_b[:,fluxEnergyChanel]],
                      'xax': xax,
                      'yax': yax,
                      'maskflux': maskflux,
                      'fontsize': 16,
                      'probe': 'A-B',
                      'level': 3,
                      'plotDir': 'rbsp/rept/',
                      'lshellOrlStar': lshellOrlStar,
                      'interpolatedFlux': False,
                      'vmin': 0.01,
                      'vmax': 5e5,
                      'figureIdentify': f'AB_HSS'
                      }
#%%


#%%
#%%

def plot_FluxvsL(**kwargs):

    for key, val in kwargs.items():
        globals()[key] = val


    out_figDir = str(Path.home().joinpath('Pictures', plotDir))

    if not os.path.exists(out_figDir) and os.path.dirname(out_figDir) != '':
        os.makedirs(out_figDir)

    xLim=[datetime.datetime.strptime(trangeXlim[0], '%Y-%m-%d'), datetime.datetime.strptime(trangeXlim[1], '%Y-%m-%d')]

    matplotlib.rc('font', size=fontsize)
    plt.close('all')
    plt.ioff()
    figprops = dict(figsize=(16,4), dpi=120)
    fig = plt.figure(**figprops)

    energChanel = f'{specEnergy[fluxEnergyChanel]}'
    if level == 2:
        om = ' -- OMNI directional'
    else:
        om = ''

    titlePanel1 = f'RBSP {probe}: ECT/REPT (L{level}) Electron Flux density, Energy {energChanel} MeV{om}'

    if figureIdentify != None:
        figureFilename = f'{figureIdentify}_RBSP{probe}_ECT-REPT_L{level}_Electron-Flux-density-Energy_{energChanel}' \
                         f'_MeV{trangeXlim[0]}_{trangeXlim[1]}.png'
    else:
        figureFilename = f'RBSP{probe}_ECT-REPT_L{level}_Electron-Flux-density-Energy_{energChanel}' \
                         f'_MeV{trangeXlim[0]}_{trangeXlim[1]}.png'


    ax = plt.axes([0.055, 0.04, 0.817, 0.88])
    plt.subplots_adjust(left=0.06, right=0.875, bottom=0.07, top=0.95)
    if interpolatedFlux == True:
        lc = ax.pcolormesh(xax, yax, maskflux, norm=colors.LogNorm(vmin=vmin, vmax = vmax), cmap='jet')
    else:
        for ff in range(len(fluxAB)):
            lc = ax.scatter(time_dt_rept[ff],l_probe[ff], 8, fluxAB[ff], norm=colors.LogNorm(vmin=vmin, vmax= vmax), cmap='jet')
    # ax.text(0.05, 0.9, '(a)', horizontalalignment='center',verticalalignment='center',
    #         fontsize=18, transform=ax.transAxes)
    ax.set_title(titlePanel1)
    ax.set_ylabel(f'{lshellOrlStar}')
    ax.set_xlabel('Time (UTC)', fontsize=15)
    ax.set_xlim(xLim[0], xLim[1])
    ax.set_xticklabels([])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    ax.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    cbar_coord = replace_at_index1(make_axes_locatable(ax).get_position(), [0,2], [0.9, 0.01])
    cbar_ax = fig.add_axes(cbar_coord)
    cbar = fig.colorbar(lc, cax=cbar_ax)
    cbar.set_label('$[cm^{-2}s^{-1}sr^{-1}MeV^{-1}]$')

    logging.info(f'saving figure at: {out_figDir}/{figureFilename}')
    plt.savefig(f'{out_figDir}/{figureFilename}')
# %%
plot_FluxvsL(**plotFluxParamsDict)
gc.collect()

#%%
exit()
ldfa = pd.DataFrame(l_rept_a, index=time_dt_rept_a, columns=['la'])
ldfb = pd.DataFrame(l_rept_b, index=time_dt_rept_b, columns=['lb'])
# %%
lldf = pd.concat([ldfa, ldfb], axis=1)
#%%

lldf['la'].values[0] is True
# %%
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, nx)

xtime, xL = np.meshgrid(np.arange(min(time), max(time), 0.01), np.arange(min(L_inerp), max(L_inerp),resolution_L))