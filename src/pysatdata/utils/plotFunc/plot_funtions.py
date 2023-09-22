import os
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import *
import datetime
from pysatdata.utils.library_functions import replace_at_index1, format_func, cutFlux_lshell
import matplotlib.colors as colors
import matplotlib.cm as cm
from pathlib import Path
from loguru import logger as logging

def plot_classicFluxSWparams(**kwargs):

    for key, val in kwargs.items():
        globals()[key] = val


    out_figDir = str(Path.home().joinpath('Pictures', plotDir))

    if not os.path.exists(out_figDir) and os.path.dirname(out_figDir) != '':
        os.makedirs(out_figDir)

    xLim=[datetime.datetime.strptime(trangeXlim[0], '%Y-%m-%d'), datetime.datetime.strptime(trangeXlim[1], '%Y-%m-%d')]

    matplotlib.rc('font', size=fontsize)
    plt.close('all')
    plt.ioff()
    figprops = dict(figsize=(16,18), dpi=120)
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


    ax = plt.axes([0.055, 0.79, 0.817, 0.18])
    plt.subplots_adjust(left=0.06, right=0.875, bottom=0.07, top=0.95)
    if interpolatedFlux == True:
        lc = ax.pcolormesh(xax, yax, maskflux, norm=colors.LogNorm(vmin=vmin, vmax = vmax), cmap='jet')
    else:
        lc = ax.scatter(time_dt_rept,l_probe, 8, flux_rept_spec[:,fluxEnergyChanel], norm=colors.LogNorm(vmin=vmin, vmax= vmax), cmap='jet')
    ax.text(0.05, 0.9, '(a)', horizontalalignment='center',verticalalignment='center',
            fontsize=18, transform=ax.transAxes)
    ax.set_title(titlePanel1)
    ax.set_ylabel(f'{lshellOrlStar}')
    # ax.set_xlabel('Time (UTC)', fontsize=15)
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

    bx = plt.axes([0.055, 0.59, 0.817, 0.18])
    cor = cm.tab10(np.linspace(0, 1, len(energyRange)))
    for il in range(len(energyRange)):
        cut_date, cutF = cutFlux_lshell(flux_rept_spec, cutLshell, il, l_probe, time_dt_rept)
        bx.semilogy(cut_date, cutF,  '-*', color=cor[il], label="{:1.1f} MeV".format(specEnergy[energyRange[il]]) )
        bx.legend(loc='upper right', bbox_to_anchor=(1.15, 1.02))
    bx.set_xlim(xLim[0], xLim[1])
    bx.text(0.05, 0.9, '(b)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=bx.transAxes)
    bx.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    bx.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    bx.get_xaxis().set_ticks_position('both')
    bx.get_yaxis().set_ticks_position('both')
    bx.set_xticklabels([])
    bx.set_title(f'Cut in flux at ${lshellOrlStar} = {cutLshell}$')
    bx.set_ylabel('$Flux / [cm^{-2}s^{-1}sr^{-1}MeV^{-1}]$')
    bx.grid()
    #
    # # Flow Speed
    cx = plt.axes([0.055, 0.48, 0.817, 0.1])
    cx.plot(time_dt_swe, flow_speed, '-', color='b')
    cx.get_xaxis().set_ticks_position('both')
    cx.get_yaxis().set_ticks_position('both')
    cx.text(0.05, 0.9, '(c)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=cx.transAxes)
    # cx.legend(loc='upper right', bbox_to_anchor=(1.12, 1.02))
    cx.set_ylabel('$V_{sw} / [Km s^{-1}]$')
    cx.set_xticklabels([])
    cx.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    cx.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    cx.set_xlim(xLim[0], xLim[1])
    cx.grid()
    #
    # ## density
    dx = plt.axes([0.055, 0.37, 0.817, 0.1])
    dx.plot(time_dt_swe, nP, '-', color='b')
    dx.get_xaxis().set_ticks_position('both')
    dx.get_yaxis().set_ticks_position('both')
    dx.text(0.05, 0.9, '(d)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=dx.transAxes)
    dx.legend(loc='upper right', bbox_to_anchor=(1.11, 1.02))
    dx.set_ylabel('$Np / [cm^{-3}]$')
    dx.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    dx.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    dx.set_xlim(xLim[0], xLim[1])
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
    ex.text(0.05, 0.9, '(e)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=ex.transAxes)
    ex.legend(loc='upper right', bbox_to_anchor=(1.14, 1.02))
    ex.set_ylabel('$IMF / [nT]$')
    ex.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    ex.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    ex.set_xlim(xLim[0], xLim[1])
    ex.set_xticklabels([])
    ex.grid()
    #
    # ## Bz
    fx = plt.axes([0.055, 0.15, 0.817, 0.1])
    fx.plot(time_dt_swe, bgse_z, '-', color='b', label="Bz_gse")
    fx.plot(time_dt_swe, b_total, '-', color='r', label="B total")
    fx.get_xaxis().set_ticks_position('both')
    fx.get_yaxis().set_ticks_position('both')
    fx.text(0.05, 0.9, '(f)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=fx.transAxes)
    fx.legend(loc='upper right', bbox_to_anchor=(1.14, 1.02))
    fx.set_ylabel('$IMF / [nT]$')
    fx.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    fx.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    fx.set_xlim(xLim[0], xLim[1])
    fx.set_xticklabels([])
    fx.grid()
    #
    #
    # ## AE index
    gx = plt.axes([0.055, 0.04, 0.817, 0.1])
    gx.plot(time_dt_swe, bgse_z, '-', color='b', label="AE Index")
    gx.get_xaxis().set_ticks_position('both')
    gx.get_yaxis().set_ticks_position('both')
    gx.text(0.05, 0.9, '(g)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=gx.transAxes)
    gx.legend(loc='upper right', bbox_to_anchor=(1.145, 1.02))
    gx.set_ylabel('$/ [nT]$')
    gx.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    gx.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    gx.set_xlim(xLim[0], xLim[1])
    gx.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    gx.grid()

    logging.info(f'saving figure at: {out_figDir}/{figureFilename}')
    plt.savefig(f'{out_figDir}/{figureFilename}')


#%%

def plot_classicSWparams(**kwargs):

    for key, val in kwargs.items():
        globals()[key] = val


    out_figDir = str(Path.home().joinpath('Pictures', plotDir))

    if not os.path.exists(out_figDir) and os.path.dirname(out_figDir) != '':
        os.makedirs(out_figDir)

    xLim=[datetime.datetime.strptime(trangeXlim[0], '%Y-%m-%d'), datetime.datetime.strptime(trangeXlim[1], '%Y-%m-%d')]

    matplotlib.rc('font', size=fontsize)
    plt.close('all')
    plt.ioff()
    figprops = dict(figsize=(16,18), dpi=120)
    fig = plt.figure(**figprops)


    if level == 2:
        om = ' -- OMNI directional'
    else:
        om = ''

    figureFilename = f'RBSP{probe}_ECT-REPT_L{level}_' \
                     f'_MeV{trangeXlim[0]}_{trangeXlim[1]}.png'

    bx = plt.axes([0.055, 0.59, 0.817, 0.1])
    bx.plot(time_dt_swe, flow_speed, '-', color='b')
    bx.set_xlim(xLim[0], xLim[1])
    bx.text(0.05, 0.9, '(b)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=bx.transAxes)
    bx.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    bx.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    bx.get_xaxis().set_ticks_position('both')
    bx.get_yaxis().set_ticks_position('both')
    bx.set_xticklabels([])
    bx.set_ylabel('$Vp [Km s^{-1}]$')
    bx.grid()
    #
    # # Flow Speed
    cx = plt.axes([0.055, 0.48, 0.817, 0.1])
    cx.plot(time_dt_swe, flow_speed, '-', color='b', label="Vsw")
    cx.get_xaxis().set_ticks_position('both')
    cx.get_yaxis().set_ticks_position('both')
    cx.text(0.05, 0.9, '(c)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=cx.transAxes)
    cx.legend(loc='upper right', bbox_to_anchor=(1.12, 1.02))
    cx.set_ylabel('$Vp [Km s^{-1}]$')
    cx.set_xticklabels([])
    cx.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    cx.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    cx.set_xlim(xLim[0], xLim[1])
    cx.grid()
    #
    # ## density
    dx = plt.axes([0.055, 0.37, 0.817, 0.1])
    dx.plot(time_dt_swe, nP, '-', color='b', label="Np")
    dx.get_xaxis().set_ticks_position('both')
    dx.get_yaxis().set_ticks_position('both')
    dx.text(0.05, 0.9, '(d)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=dx.transAxes)
    dx.legend(loc='upper right', bbox_to_anchor=(1.11, 1.02))
    dx.set_ylabel('$\\# / [cm^{-3}]$')
    dx.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    dx.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    dx.set_xlim(xLim[0], xLim[1])
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
    ex.text(0.05, 0.9, '(e)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=ex.transAxes)
    ex.legend(loc='upper right', bbox_to_anchor=(1.14, 1.02))
    ex.set_ylabel('$\\# / [nT]$')
    ex.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    ex.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    ex.set_xlim(xLim[0], xLim[1])
    ex.set_xticklabels([])
    ex.grid()
    #
    # ## Bz
    fx = plt.axes([0.055, 0.15, 0.817, 0.1])
    fx.plot(time_dt_swe, bgse_z, '-', color='b', label="Bz_gse")
    fx.plot(time_dt_swe, b_total, '-', color='r', label="B total")
    fx.get_xaxis().set_ticks_position('both')
    fx.get_yaxis().set_ticks_position('both')
    fx.text(0.05, 0.9, '(f)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=fx.transAxes)
    fx.legend(loc='upper right', bbox_to_anchor=(1.14, 1.02))
    fx.set_ylabel('$\\# / [nT]$')
    fx.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    fx.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    fx.set_xlim(xLim[0], xLim[1])
    fx.set_xticklabels([])
    fx.grid()
    #
    #
    # ## AE index
    gx = plt.axes([0.055, 0.04, 0.817, 0.1])
    gx.plot(time_dt_swe, ae_index, '-', color='b', label="AE Index")
    gx.get_xaxis().set_ticks_position('both')
    gx.get_yaxis().set_ticks_position('both')
    gx.text(0.05, 0.9, '(g)', horizontalalignment='center', verticalalignment='center',
            fontsize=18, transform=gx.transAxes)
    gx.legend(loc='upper right', bbox_to_anchor=(1.145, 1.02))
    gx.set_ylabel('$\\# / [nT]$')
    gx.tick_params(direction='in', length=10, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='major')
    gx.tick_params(direction='in', length=7, width=0.7, colors='k',
                   grid_color='k', grid_alpha=0.5, which='minor')
    gx.set_xlim(xLim[0], xLim[1])
    gx.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    gx.grid()

    logging.info(f'saving figure at: {out_figDir}/{figureFilename}')
    plt.savefig(f'{out_figDir}/{figureFilename}')


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


    ax = plt.axes([0.055, 0.79, 0.817, 0.18])
    plt.subplots_adjust(left=0.06, right=0.875, bottom=0.07, top=0.95)
    if interpolatedFlux == True:
        lc = ax.pcolormesh(xax, yax, maskflux, norm=colors.LogNorm(vmin=vmin, vmax = vmax), cmap='jet')
    else:
        for ff in fluxAB:
            lc = ax.scatter(time_dt_rept,l_probe, 8, ff, norm=colors.LogNorm(vmin=vmin, vmax= vmax), cmap='jet')
    ax.text(0.05, 0.9, '(a)', horizontalalignment='center',verticalalignment='center',
            fontsize=18, transform=ax.transAxes)
    ax.set_title(titlePanel1)
    ax.set_ylabel(f'{lshellOrlStar}')
    # ax.set_xlabel('Time (UTC)', fontsize=15)
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