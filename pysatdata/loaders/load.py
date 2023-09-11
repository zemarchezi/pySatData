from pysatdata.utils.dailynames import dailynames
from pysatdata.utils.library_functions import testRemoteDir, testFiles
from pysatdata.utils.download import download
from pysatdata.resources.config import openConfigFile
from pysatdata.satellites.read_goes import *
from pysatdata.satellites.read_rbsp import *
from pysatdata.satellites.read_ace import *
from pysatdata.satellites.read_ace_merged import *
from pysatdata.satellites.read_themis import *
from pathlib import Path
from loguru import logger as logging
import os
# %%

def load_sat(trange: list=['2013-11-5', '2013-11-6'],
             satellite: str='goes',
             probe: list=[],
             level: str='2',
             instrument: str='magn',
             datatype: str='sectors',
             suffix: str='',
             cadence: str='4sec', # for EMFISIS and ACE mag data
             coord: str='sm', # for EMFISIS mag data
             rel: str='rel04', # for ECT data
             time: str='time', # For GOES satellites
             get_support_data: bool = False,
             varformat=None,
             varnames: list=[],
             downloadonly: bool=False,
             notplot: bool=False,
             no_update: bool=False,
             time_clip: bool=False,
             usePandas: bool = True,
             usePyTplot: bool = False,
             testRemotePath: bool=False,
             searchFilesFirst: bool=False,
             config_file: str = './pysatdata/resources/config_file.json'):


    global remote_path, out_files, pathformat, tvars

    # config_file = openConfigFile(path=config_file).config_file_sat
    config_file = openConfigFile().config_file_sat
    # set the download data directories
    if os.environ.get('DATA_DIR'):
        config_file[satellite]['local_data_dir'] = os.sep.join([os.environ['DATA_DIR'], f'{satellite}'])

    if os.environ.get(f'{satellite.upper()}_DATA_DIR'):
        config_file[satellite]['local_data_dir'] = os.environ[f'{satellite.upper()}_DATA_DIR']

    # local data path
    if config_file[satellite]['local_data_dir'] != "sat_data/":
        local_path = f"/{config_file[satellite]['local_data_dir']}{ satellite}"
    else:
        local_path = str(Path.home().joinpath(config_file[satellite]['local_data_dir'], satellite))
    logging.info(f'Local Download Path: {local_path}')



    for prb in probe:
        logging.warning("Selecting the sub path key")
        # test remote Path
        if testRemotePath:
            remote_path, subpathKey, filenameKey, datatype = \
                testRemoteDir(config_file, satellite, prb, instrument, level, datatype)
        else:
            remote_path = config_file[satellite]['remote_data_dir']
            subpathKey = 'subpath'
            filenameKey = 'filename'
        logging.info(f'Remotepath: {remote_path}')
        # subpathKey = "altern_subpath"
        # filenameKey = "altern_filename"
        if instrument == "MagEphem":
            remote_path = config_file[satellite]['remote_subpath'][str(prb)][instrument]["secondRemoteDir"]
        # # datatype = 'pitchangle'


        subpathformat = config_file[satellite]['remote_subpath'][str(prb)][instrument][datatype][subpathKey]
        subpathformat = eval(f"f'{subpathformat}'")
        sat_filename = config_file[satellite]['remote_subpath'][str(prb)][instrument][datatype][filenameKey]
        sat_filename = eval(f"f'{sat_filename}'")
        pathformat = f"{subpathformat}{sat_filename}"
        logging.info(f'Remote_file_path: {pathformat}')

        # find the full remote path names using the trange
        remote_names = dailynames(file_format=pathformat, trange=trange)
        out_files = []
        if searchFilesFirst:
            remote_names, files = testFiles(local_path, remote_names)
            if files is not None:
                for file in files:
                    out_files.append(file)
        # download the files
        files = download(remote_file=remote_names, remote_path=remote_path, local_path=local_path, no_download=no_update)
        if files is not None:
            for file in files:
                out_files.append(file)

    out_files = sorted(out_files)

    if downloadonly:
        return out_files

    if satellite in ['goes16_17', 'goes13_15']:

        tvars = readData_goes(out_files, usePyTplot, usePandas, suffix, time)

    if satellite in ['rbsp', 'erg']:
        tvars = readData_rbsp(out_files, usePyTplot, usePandas, suffix, get_support_data,
                              varformat, varnames, notplot)
    if satellite in ['ace', 'omni']:
        if instrument in ['merged']:
            tvars = readData_ace_merged(out_files, usePyTplot, usePandas, suffix, get_support_data,
                              varformat, varnames, notplot)
        else:
            tvars = readData_ace(out_files, usePyTplot, usePandas, suffix, get_support_data,
                              varformat, varnames, notplot)
    if satellite == 'themis':
        tvars = readData_themis(out_files, usePyTplot, usePandas, suffix, get_support_data,
                              varformat, varnames, notplot)
    if time_clip:
        for new_var in tvars:
            tclip(new_var, trange[0], trange[1], suffix='')

    return tvars

# if satellite == 'ace':
#     remote_path = config_file[satellite]['remote_data_dir']
#     logging.info(f'Remotepath: {remote_path}')
#
#     if instrument == 'fgm':
#         pathformat = f'mag/level_2_cdaweb/mfi_h{level}/%Y/ac_{datatype}_mfi_%Y%m%d_v??.cdf'
#     elif instrument == 'swe':
#         pathformat = 'swepam/level_2_cdaweb/swe_' + datatype + '/%Y/ac_' + datatype + '_swe_%Y%m%d_v??.cdf'
#     elif instrument == 'epm':
#         pathformat = 'epam/level_2_cdaweb/epm_' + datatype + '/%Y/ac_' + datatype + '_epm_%Y%m%d_v??.cdf'
#     elif instrument == 'cris':
#         pathformat = 'cris/level_2_cdaweb/cris_' + datatype + '/%Y/ac_' + datatype + '_cris_%Y%m%d_v??.cdf'
#     elif instrument == 'sis':
#         pathformat = 'sis/level_2_cdaweb/sis_' + datatype + '/%Y/ac_' + datatype + '_sis_%Y%m%d_v??.cdf'
#     elif instrument == 'ule':
#         pathformat = 'uleis/level_2_cdaweb/ule_' + datatype + '/%Y/ac_' + datatype + '_ule_%Y%m%d_v??.cdf'
#     elif instrument == 'sep':
#         pathformat = 'sepica/level_2_cdaweb/sep_' + datatype + '/%Y/ac_' + datatype + '_sep_%Y%m%d_v??.cdf'
#     elif instrument == 'swics':
#         filename_dtype = datatype.split('_')[1] + '_' + datatype.split('_')[0]
#         pathformat = 'swics/level_2_cdaweb/' + datatype + '/%Y/ac_' + filename_dtype + '_%Y%m%d_v??.cdf'
#
# "remote_subpath": {
#       "ace": {
#         "swepam": {
#           "level_2_cdaweb": {
#             "subpath": "rbsp{prb}/{level}/{instrument}/{datatype}/%Y/",
#             "filename": "rbsp-{prb}_{datatype}_{instrument}-{level}_%Y%m%d_v*.cdf",
#           },
#         },
#       },
#     }
