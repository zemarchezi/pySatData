from pysatdata.utils.dailynames import dailynames
from pysatdata.utils.download import download
from pysatdata.satellites.read_goes import *
from pysatdata.satellites.read_rbsp import *
from pathlib import Path
import glob
import logging

# %%

def load_sat(trange: list=['2013-11-5', '2013-11-6'],
             satellite: str='goes',
             probe: list=[],
             level: str='2',
             instrument: str='magn',
             datatype: str='hires',
             suffix: str='',
             cadence='4sec', # for EMFISIS mag data
             coord='sm', # for EMFISIS mag data
             rel='rel04', # for ECT data
             get_support_data: bool = False,
             varformat=None,
             varnames: list=[],
             downloadonly: bool=False,
             notplot: bool=False,
             no_update: bool=False,
             time_clip: bool=False,
             usePandas: bool = True,
             usePyTplot: bool = False,
             config_file: dict = {}):


    global remote_path, out_files, pathformat, tvars

    local_path = str(Path.home().joinpath(config_file[satellite]['local_data_dir'], satellite))
    logging.info(f'Local Download Path: {local_path}')

    for prb in probe:
        if satellite == 'goes':
            if int(prb) > 15:
                remote_path = config_file[satellite]['remote_data_dir']
                logging.info(f'Remotepath: {remote_path}')
                pathformat = f"goes{prb}/l{level}/data/{instrument}-l{level}-{datatype}/%Y/%m/dn_{instrument}-l{level}-{datatype}_g{prb}_d%Y%m%d_*.nc"
            else:
                fullavgpath = ['full', 'avg']
                goes_path_dir = fullavgpath[datatype == '1min' or datatype == '5min']
                remote_path = goes_path_dir + '/%Y/%m/goes' + str(prb) + '/netcdf/'
                if instrument == 'fgm':
                    if datatype == '512ms': # full, unaveraged data
                        pathformat = remote_path + 'g' + str(prb) + '_magneto_512ms_%Y%m%d_%Y%m%d.nc'
                    elif datatype == '1min': # 1 min averages
                        pathformat = remote_path + 'g' + str(prb) + '_magneto_1m_%Y%m01_%Y%m??.nc'
                    elif datatype == '5min': # 5 min averages
                        pathformat = remote_path + 'g' + str(prb) + '_magneto_5m_%Y%m01_%Y%m??.nc'

        if satellite == 'rbsp':
            remote_path = config_file[satellite]['remote_data_dir']
            logging.info(f'Remotepath: {remote_path}')
            if instrument == 'emfisis':
                subpathformat = f"rbsp{prb}/{level}/{instrument}/{datatype}/{cadence}/{coord}/%Y/"
                pathformat = f"{subpathformat}rbsp-{prb}_{datatype}_{cadence}-{coord}_{instrument}-{level}_%Y%m%d_v*.cdf"
            elif instrument == 'rbspice':
                subpathformat = f"rbsp{prb}/{level}/{instrument}/{datatype}/%Y/"
                pathformat = f"{subpathformat}rbsp-{prb}-{instrument}_lev-{str(level)[-1]}_{datatype}_%Y%m%d_v*.cdf"
            elif instrument == 'efw':
                subpathformat = f"rbsp{prb}/{level}/{instrument}/%Y/"
                pathformat = f"{subpathformat}rbsp{prb}_{instrument}-{level}_%Y%m%d_v??.cdf"
            elif instrument == 'mageis':
                subpathformat = f"rbsp{prb}/{level}/ect/{instrument}/sectors/{rel}/%Y/"
                pathformat = f"{subpathformat}rbsp{prb}_{rel}_ect-mageis-{level}_%Y%m%d_v*.cdf"
            elif instrument == 'hope':
                subpathformat = f"rbsp{prb}/{level}/ect/{instrument}/{datatype}/{rel}/%Y/"
                if datatype == 'moments':
                    pathformat = f"{subpathformat}rbsp{prb}_{rel}_ect-hope-mom-{level}_%Y%m%d_v*.cdf"
                elif datatype == 'pitchangle':
                    pathformat = f"{subpathformat}rbsp{prb}_{rel}_ect-hope-pa-{level}_%Y%m%d_v*.cdf"
                elif datatype == 'spinaverage':
                    pathformat = f"{subpathformat}rbsp{prb}_{rel}_ect-hope-sci-{level}sa_%Y%m%d_v*.cdf"
                elif datatype == 'sectors':
                    pathformat = f"{subpathformat}rbsp{prb}_{rel}_ect-hope-sci-{level}_%Y%m%d_v*.cdf"
            elif instrument == 'rept':
                subpathformat = f"rbsp{prb}/{level}/ect/{instrument}/sectors/{rel}/%Y/"
                pathformat = f"{subpathformat}rbsp{prb}_{rel}_ect-rept-sci-{level}_%Y%m%d_v*.cdf"
            elif instrument == 'rps':
                subpathformat = f"rbsp{prb}/{level}/psbr/{datatype}/%Y/"
                if datatype == 'rps-1min':
                    pathformat = f"{subpathformat}rbsp{prb}_{level}-1min_psbr-rps_%Y%m%d_v*.cdf"
                elif datatype == 'rps':
                    pathformat = f"{subpathformat}rbsp{prb}_{level}_psbr-rps_%Y%m%d_v*.cdf"
        # find the full remote path names using the trange
        remote_names = dailynames(file_format=pathformat, trange=trange)

        out_files = []

        files = download(remote_file=remote_names, remote_path=remote_path, local_path=local_path, no_download=no_update)
        if files is not None:
            for file in files:
                out_files.append(file)

    out_files = sorted(out_files)

    if downloadonly:
        return out_files

    if satellite == 'goes':
        tvars = readData_goes(out_files, usePyTplot, usePandas, suffix, time='time')

    if satellite == 'rbsp':
        tvars = readData_rbsp(out_files, usePyTplot, usePandas, suffix, get_support_data,
                              varformat, varnames, notplot)
    if time_clip:
        for new_var in tvars:
            tclip(new_var, trange[0], trange[1], suffix='')

    return tvars
