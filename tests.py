from pysatdata.loaders.load import *
import json
import pytplot
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import *
import datetime
import pytz
from pysatdata.utils.library_functions import replace_at_index1, format_func, cutFlux_lshell
from pysatdata.utils.interpolate_flux_rbsp import *
#%%

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

xax, yax, maskflux = interpolateFluxRbsp(flux_rept_spec[:,1], l_rept, time_dt_rept)

#%%
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
matplotlib.rc('font', size=15)
plt.close('all')
plt.ioff()
figprops = dict(figsize=(16,18), dpi=120)
fig = plt.figure(**figprops)

ax = plt.axes([0.055, 0.79, 0.817, 0.18])
plt.subplots_adjust(left=0.06, right=0.875, bottom=0.07, top=0.95)
lc = ax.pcolormesh(xax, yax, maskflux, norm=colors.LogNorm(vmin=0.2, vmax = 5e5), cmap='viridis')
ax.set_title('RBSP A: ECT/REPT (L2) Electron Flux density, Energy 2.10 MeV -- OMNI directional')
ax.set_ylabel('L-Shell', fontsize=13)
# ax.set_xlabel('Time (UTC)', fontsize=15)
ax.set_xlim(datetime.datetime.strptime(trange0[0], '%Y-%m-%d'), datetime.datetime.strptime(trange0[1], '%Y-%m-%d'))
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
# ax[0].xaxis.set_major_formatter(plt.FuncFormatter(format_func))
# ax.grid()

bx = plt.axes([0.055, 0.59, 0.817, 0.18])
cor = cm.tab10(np.linspace(0, 1, 4))
for il in range(0,4):
    cut_date, cutF = cutFlux_lshell(flux_rept_spec, 5, il, l_rept, time_dt_rept)
    bx.semilogy(cut_date, cutF,  '-*', color=cor[il], label="{:1.1f} MeV".format(spec[il]) )
    bx.legend(loc='upper right', bbox_to_anchor=(1.15, 1.02))
bx.set_xlim(datetime.datetime.strptime(trange0[0], '%Y-%m-%d'), datetime.datetime.strptime(trange0[1], '%Y-%m-%d'))
bx.tick_params(direction='in', length=10, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='major')
bx.tick_params(direction='in', length=7, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='minor')
bx.get_xaxis().set_ticks_position('both')
bx.get_yaxis().set_ticks_position('both')
bx.set_xticklabels([])
bx.set_title('Cut in flux at $L = 5$')
bx.set_ylabel('$\\# / [cm^{-2}s^{-1}sr^{-1}MeV^{-1}]$')
bx.grid()
#
# # Flow Speed
cx = plt.axes([0.055, 0.48, 0.817, 0.1])
cx.plot(time_dt_swe, flow_speed, '-', color='b', label="Vsw")
cx.get_xaxis().set_ticks_position('both')
cx.get_yaxis().set_ticks_position('both')
cx.legend(loc='upper right', bbox_to_anchor=(1.12, 1.02))
cx.set_ylabel('$\\# / [Km s^{-1}]$')
cx.set_xticklabels([])
cx.tick_params(direction='in', length=10, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='major')
cx.tick_params(direction='in', length=7, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='minor')
cx.set_xlim(datetime.datetime.strptime(trange0[0], '%Y-%m-%d'), datetime.datetime.strptime(trange0[1], '%Y-%m-%d'))
cx.grid()
#
# ## density
dx = plt.axes([0.055, 0.37, 0.817, 0.1])
dx.plot(time_dt_swe, nP, '-', color='b', label="Np")
dx.get_xaxis().set_ticks_position('both')
dx.get_yaxis().set_ticks_position('both')
dx.legend(loc='upper right', bbox_to_anchor=(1.11, 1.02))
dx.set_ylabel('$\\# / [cm^{-3}]$')
dx.tick_params(direction='in', length=10, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='major')
dx.tick_params(direction='in', length=7, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='minor')
dx.set_xlim(datetime.datetime.strptime(trange0[0], '%Y-%m-%d'), datetime.datetime.strptime(trange0[1], '%Y-%m-%d'))
dx.set_xticklabels([])
dx.grid()
#
#
# ## imf
ex = plt.axes([0.055, 0.26, 0.817, 0.1])
ex.plot(time_dt_swe, bgse_x, '-', color='b', label="Bx_gse")
ex.plot(time_dt_swe, bgse_y, '-', color='r', label="By_gse")
ex.get_xaxis().set_ticks_position('both')
ex.get_yaxis().set_ticks_position('both')
ex.legend(loc='upper right', bbox_to_anchor=(1.14, 1.02))
ex.set_ylabel('$\\# / [nT]$')
ex.tick_params(direction='in', length=10, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='major')
ex.tick_params(direction='in', length=7, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='minor')
ex.set_xlim(datetime.datetime.strptime(trange0[0], '%Y-%m-%d'), datetime.datetime.strptime(trange0[1], '%Y-%m-%d'))
ex.set_xticklabels([])
ex.grid()
#
# ## Bz
fx = plt.axes([0.055, 0.15, 0.817, 0.1])
fx.plot(time_dt_swe, bgse_z, '-', color='b', label="Bz_gse")
fx.plot(time_dt_swe, b_total, '-', color='r', label="B total")
fx.get_xaxis().set_ticks_position('both')
fx.get_yaxis().set_ticks_position('both')
fx.legend(loc='upper right', bbox_to_anchor=(1.14, 1.02))
fx.set_ylabel('$\\# / [nT]$')
fx.tick_params(direction='in', length=10, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='major')
fx.tick_params(direction='in', length=7, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='minor')
fx.set_xlim(datetime.datetime.strptime(trange0[0], '%Y-%m-%d'), datetime.datetime.strptime(trange0[1], '%Y-%m-%d'))
fx.set_xticklabels([])
fx.grid()
#
#
# ## AE index
gx = plt.axes([0.055, 0.04, 0.817, 0.1])
gx.plot(time_dt_swe, ae_index, '-', color='b', label="AE Index")
gx.get_xaxis().set_ticks_position('both')
gx.get_yaxis().set_ticks_position('both')
gx.legend(loc='upper right', bbox_to_anchor=(1.145, 1.02))
gx.set_ylabel('$\\# / [nT]$')
gx.tick_params(direction='in', length=10, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='major')
gx.tick_params(direction='in', length=7, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='minor')
gx.set_xlim(datetime.datetime.strptime(trange0[0], '%Y-%m-%d'), datetime.datetime.strptime(trange0[1], '%Y-%m-%d'))
gx.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
gx.grid()
plt.savefig('testefluxLcut.png')

#%%
#
# cut_date, cutF = cutFlux_lshell(flux_rept_spec,5, 1, l_rept, time_dt_rept)
#
#
# #%%
# fig, ax = plt.subplots(1, 1, figsize=(16,4), dpi=120)
# cor = cm.tab10(np.linspace(0, 1, 4))
# for il in range(0,4):
#     cut_date, cutF = cutFlux_lshell(flux_rept_spec, 5, il, l_rept, time_dt_rept)
#     ax.semilogy(cut_date, cutF,  '-*', color=cor[il], label="{:1.1f} MeV".format(spec[il]) )
#     plt.legend()
# ax.set_xlim(datetime.datetime.strptime(trange0[0], '%Y-%m-%d'), datetime.datetime.strptime(trange0[1], '%Y-%m-%d'))
# ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
# plt.grid()
# plt.savefig('cutTest.png')
#
#
# #%%
#
# date_range = trange
# fig, ax = plt.subplots(3, 1, figsize=(16,12))
#
# # matplotlib.rc('xtick', labelsize=15)
# # matplotlib.rc('ytick', labelsize=15)
#
# lc = ax[0].scatter(time_dt_rept, l_rept, 7, flux_rept_spec[:,1], norm=colors.LogNorm(vmin=0.5, vmax = 5e6), cmap='jet')
# ax[0].set_title('RBSP A: ECT/REPT (L3) Intensity: Electrons Pitch Angle 90.0-100.6 Deg 2.10 MeV')
# ax[0].set_ylabel('L-Shell', fontsize=13)
# # ax.set_xlabel('Time (UTC)', fontsize=15)
# ax[0].set_xlim(np.min(time_dt_rept), np.max(time_dt_rept))
# ax[0].set_xticklabels([])
# cbar_coord = replace_at_index1(make_axes_locatable(ax[0]).get_position(), [0,2], [0.92, 0.01])
# cbar_ax = fig.add_axes(cbar_coord)
# cbar = fig.colorbar(lc, cax=cbar_ax)
# cbar.set_label('$[cm^{-2}s^{-1}sr^{-1}MeV^{-1}]$')
# ax[0].grid()
#
# lc2 = ax[1].scatter(time_dt_rept, l_rept, 7, flux_rept_spec[:,1], norm=colors.LogNorm(vmin=1e3, vmax = 1e6), cmap='jet')
# ax[1].set_title('RBSP A: ECT/MAGEIS (L3) Intensity: Electrons Pitch Angle 73.64-90.00 Deg 33 keV')
# ax[1].set_ylabel('L-Shell', fontsize=13)
# # ax.set_xlabel('Time (UTC)', fontsize=15)
# ax[1].set_xlim(np.min(time_dt_rept), np.max(time_dt_rept))
# ax[1].set_xticklabels([])
# cbar_coord = replace_at_index1(make_axes_locatable(ax[1]).get_position(), [0,2], [0.92, 0.01])
# cbar_ax = fig.add_axes(cbar_coord)
# cbar = fig.colorbar(lc2, cax=cbar_ax)
# cbar.set_label('$[cm^{-2}s^{-1}sr^{-1}keV^{-1}]$')
# ax[1].grid()
#
#
# lc3 = ax[2].scatter(time_dt_rept, l_rept, 7, flux_rept_spec[:,1], norm=colors.LogNorm(vmin=1e5, vmax = 1e10), cmap='jet')
# ax[2].set_title('RBSP A: ECT/HOPE (L3) Intensity: Electrons Pitch Angle 90.0 Deg 9259.82 eV')
# ax[2].set_ylabel('L-Shell', fontsize=13)
# # ax.set_xlabel('Time (UTC)', fontsize=15)
# ax[2].set_xlim(np.min(time_dt_rept), np.max(time_dt_rept))
# ax[2].xaxis.set_major_formatter(plt.FuncFormatter(format_func))
# cbar_coord = replace_at_index1(make_axes_locatable(ax[2]).get_position(), [0,2], [0.92, 0.01])
# cbar_ax = fig.add_axes(cbar_coord)
# cbar = fig.colorbar(lc3, cax=cbar_ax)
# cbar.set_label('$[cm^{-2}s^{-1}sr^{-1}eV^{-1}]$')
# ax[2].grid()
#
# plt.savefig(f"Flux_{date_range[0]}_to_{date_range[1]}.png")
# # plt.show()
#
#
# #%%
#
#
# # matplotlib.rc('xtick', labelsize=15)
# # matplotlib.rc('ytick', labelsize=15)
