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

with open('/home/jose/python_projects/pySatData/pysatdata/utils/resources/config_file.json', 'r') as f:
    config_file_sat = json.load(f)


trange0=['2016-03-06', '2016-03-10']
trange=['2016-03-05', '2016-03-11']
varss_rept = load_sat(trange=trange, satellite='rbsp',
                     probe=['a'], level='l2', rel='rel03',
                     instrument='rept',
                     config_file=config_file_sat, downloadonly=False,
                     usePandas=False, usePyTplot=True)
#%%
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

#%%
xax, yax, maskflux = interpolateFluxRbsp(flux_rept_spec[:,1], l_rept, time_dt_rept)

#%%


fig, ax = plt.subplots(2, 1, figsize=(16,8), dpi=120)
plt.subplots_adjust(left=0.06, bottom=0.07, top=0.95)
lc = ax[0].pcolormesh(xax, yax, maskflux, norm=colors.LogNorm(vmin=0.2, vmax = 5e5), cmap='viridis')
ax[0].set_title('RBSP A: ECT/REPT (L2) Electron Flux density, Energy 2.10 MeV -- OMNI directional')
ax[0].set_ylabel('L-Shell', fontsize=13)
# ax.set_xlabel('Time (UTC)', fontsize=15)
ax[0].set_xlim(datetime.datetime.strptime(trange0[0], '%Y-%m-%d'), datetime.datetime.strptime(trange0[1], '%Y-%m-%d'))
ax[0].set_xticklabels([])
ax[0].xaxis.set_minor_locator(AutoMinorLocator())
ax[0].tick_params(direction='in', length=10, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='major')
ax[0].tick_params(direction='in', length=7, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='minor')
cbar_coord = replace_at_index1(make_axes_locatable(ax[0]).get_position(), [0,2], [0.92, 0.01])
cbar_ax = fig.add_axes(cbar_coord)
cbar = fig.colorbar(lc, cax=cbar_ax)
cbar.set_label('$[cm^{-2}s^{-1}sr^{-1}MeV^{-1}]$')
# ax[0].xaxis.set_major_formatter(plt.FuncFormatter(format_func))
# ax.grid()


cor = cm.tab10(np.linspace(0, 1, 4))
for il in range(0,4):
    cut_date, cutF = cutFlux_lshell(flux_rept_spec, 5, il, l_rept, time_dt_rept)
    ax[1].semilogy(cut_date, cutF,  '-*', color=cor[il], label="{:1.1f} MeV".format(spec[il]) )
    ax[1].legend(loc='upper right', bbox_to_anchor=(1.1, 1.02))
ax[1].set_xlim(datetime.datetime.strptime(trange0[0], '%Y-%m-%d'), datetime.datetime.strptime(trange0[1], '%Y-%m-%d'))
ax[1].tick_params(direction='in', length=10, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='major')
ax[1].tick_params(direction='in', length=7, width=0.7, colors='k',
               grid_color='k', grid_alpha=0.5, which='minor')
ax[1].get_xaxis().set_ticks_position('both')
ax[1].get_yaxis().set_ticks_position('both')
ax[1].set_title('Cut in flux at $L = 5$')
ax[1].set_ylabel('Flux $[cm^{-2}s^{-1}sr^{-1}MeV^{-1}]$')
ax[1].xaxis.set_major_formatter(plt.FuncFormatter(format_func))
ax[1].grid()


plt.savefig('testefluxLcut.png')

#%%

cut_date, cutF = cutFlux_lshell(flux_rept_spec,5, 1, l_rept, time_dt_rept)


#%%
fig, ax = plt.subplots(1, 1, figsize=(16,4), dpi=120)
cor = cm.tab10(np.linspace(0, 1, 4))
for il in range(0,4):
    cut_date, cutF = cutFlux_lshell(flux_rept_spec, 5, il, l_rept, time_dt_rept)
    ax.semilogy(cut_date, cutF,  '-*', color=cor[il], label="{:1.1f} MeV".format(spec[il]) )
    plt.legend()
ax.set_xlim(datetime.datetime.strptime(trange0[0], '%Y-%m-%d'), datetime.datetime.strptime(trange0[1], '%Y-%m-%d'))
ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
plt.grid()
plt.savefig('cutTest.png')


#%%

date_range = trange
fig, ax = plt.subplots(3, 1, figsize=(16,12))

# matplotlib.rc('xtick', labelsize=15)
# matplotlib.rc('ytick', labelsize=15)

lc = ax[0].scatter(time_dt_rept, l_rept, 7, flux_rept_spec[:,1], norm=colors.LogNorm(vmin=0.5, vmax = 5e6), cmap='jet')
ax[0].set_title('RBSP A: ECT/REPT (L3) Intensity: Electrons Pitch Angle 90.0-100.6 Deg 2.10 MeV')
ax[0].set_ylabel('L-Shell', fontsize=13)
# ax.set_xlabel('Time (UTC)', fontsize=15)
ax[0].set_xlim(np.min(time_dt_rept), np.max(time_dt_rept))
ax[0].set_xticklabels([])
cbar_coord = replace_at_index1(make_axes_locatable(ax[0]).get_position(), [0,2], [0.92, 0.01])
cbar_ax = fig.add_axes(cbar_coord)
cbar = fig.colorbar(lc, cax=cbar_ax)
cbar.set_label('$[cm^{-2}s^{-1}sr^{-1}MeV^{-1}]$')
ax[0].grid()

lc2 = ax[1].scatter(time_dt_rept, l_rept, 7, flux_rept_spec[:,1], norm=colors.LogNorm(vmin=1e3, vmax = 1e6), cmap='jet')
ax[1].set_title('RBSP A: ECT/MAGEIS (L3) Intensity: Electrons Pitch Angle 73.64-90.00 Deg 33 keV')
ax[1].set_ylabel('L-Shell', fontsize=13)
# ax.set_xlabel('Time (UTC)', fontsize=15)
ax[1].set_xlim(np.min(time_dt_rept), np.max(time_dt_rept))
ax[1].set_xticklabels([])
cbar_coord = replace_at_index1(make_axes_locatable(ax[1]).get_position(), [0,2], [0.92, 0.01])
cbar_ax = fig.add_axes(cbar_coord)
cbar = fig.colorbar(lc2, cax=cbar_ax)
cbar.set_label('$[cm^{-2}s^{-1}sr^{-1}keV^{-1}]$')
ax[1].grid()


lc3 = ax[2].scatter(time_dt_rept, l_rept, 7, flux_rept_spec[:,1], norm=colors.LogNorm(vmin=1e5, vmax = 1e10), cmap='jet')
ax[2].set_title('RBSP A: ECT/HOPE (L3) Intensity: Electrons Pitch Angle 90.0 Deg 9259.82 eV')
ax[2].set_ylabel('L-Shell', fontsize=13)
# ax.set_xlabel('Time (UTC)', fontsize=15)
ax[2].set_xlim(np.min(time_dt_rept), np.max(time_dt_rept))
ax[2].xaxis.set_major_formatter(plt.FuncFormatter(format_func))
cbar_coord = replace_at_index1(make_axes_locatable(ax[2]).get_position(), [0,2], [0.92, 0.01])
cbar_ax = fig.add_axes(cbar_coord)
cbar = fig.colorbar(lc3, cax=cbar_ax)
cbar.set_label('$[cm^{-2}s^{-1}sr^{-1}eV^{-1}]$')
ax[2].grid()

plt.savefig(f"Flux_{date_range[0]}_to_{date_range[1]}.png")
# plt.show()


#%%


# matplotlib.rc('xtick', labelsize=15)
# matplotlib.rc('ytick', labelsize=15)
