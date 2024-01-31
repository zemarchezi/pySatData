#%%
from pysatdata.loaders.load import *
from pysatdata.utils.interpolate_flux_rbsp import * 
from pysatdata.utils.plotFunc.plot_funtions import *
from pysatdata.utils.library_functions import *
#%%
trange= ['2012-11-01','2012-11-03']
downloadonly = True
paramLoadrept = {"satellite": 'rbsp', "probe": ['a'], "level": '3', "rel": "rel03",
                    "instrument": 'rept', "datatype": "sectors"}
#%%
varss_rept = load_sat(trange=trange, satellite=paramLoadrept['satellite'],
                        probe=paramLoadrept['probe'], level='2', 
                        rel=paramLoadrept['rel'], instrument=paramLoadrept['instrument'],
                        datatype=paramLoadrept['datatype'],testRemotePath=False,
                        downloadonly=downloadonly, 
                        searchFilesFirst=True, usePandas=False, usePyTplot=True)
#%%

trange_plot=['2012-07-03', '2012-07-05'] # time range for ploting
trange=['2012-07-02', '2012-07-06'] 
paramLoadSat = {"satellite": 'dmsp', "probe": 'f17',
                "instrument": 'ssies', "datatype": 'thermal-plasma'}

varss_dmsp = load_sat(trange=trange, satellite=paramLoadSat["satellite"],
                     probe=[paramLoadSat["probe"]],
                     instrument=paramLoadSat["instrument"],datatype=paramLoadSat["datatype"], downloadonly=False,
                     testRemotePath=True,
                     usePandas=False, usePyTplot=True)
print(varss_dmsp)
# for n, dd in enumerate(dataT.index):


def downloadData(trange, config_file):


    config_file_sat = config_file



    lshellOrlStar = 'L-Star'
    probe = 'a'

    paramLoadEfw = {"satellite": 'rbsp', "probe": ['a', 'b'], "level": '2', "rel": "rel03",
                    "instrument": 'efw', "datatype": 'esvy_despun',
                    "varnames" : ['efield_mgse', 'lshell']}
    paramLoadEm = {"satellite": 'rbsp', "probe": ['a', 'b'], "level": '3', "rel": "rel03",
                    "instrument": 'emfisis', "datatype": "magnetometer", "coord": "gse",
                    "cadence": "1sec", "varnames": []}
    paramLoadrept = {"satellite": 'rbsp', "probe": ['a', 'b'], "level": '3', "rel": "rel03",
                    "instrument": 'rept', "datatype": "sectors"}
    paramLoadmageis = {"satellite": 'rbsp', "probe": ['a', 'b'], "level": '3', "rel": "rel04",
                    "instrument": 'rept', "datatype": "sectors"}
    #%%

    downloadonly = True

    varss_emfisis = load_sat(trange=trange, satellite=paramLoadEm['satellite'],
                            probe=[probe], rel='rel03', level=paramLoadEm['level'],
                            instrument=paramLoadEm['instrument'], datatype=paramLoadEm['datatype'],
                            cadence=paramLoadEm['cadence'], coord=paramLoadEm['coord'],
                            varnames=paramLoadEm['varnames'], downloadonly=downloadonly,
                            usePandas=True, usePyTplot=False)



    pytplot.del_data()
    varss_rept = load_sat(trange=trange, satellite=paramLoadrept['satellite'],
                        probe=paramLoadrept['probe'], level=paramLoadrept['level'], 
                        rel=paramLoadrept['rel'], instrument=paramLoadrept['instrument'],
                        datatype=paramLoadrept['datatype'],
                        downloadonly=downloadonly, 
                        usePandas=False, usePyTplot=True)
    pytplot.del_data()
    varss_rept = load_sat(trange=trange, satellite=paramLoadrept['satellite'],
                        probe=paramLoadrept['probe'], level='2', 
                        rel=paramLoadrept['rel'], instrument=paramLoadrept['instrument'],
                        datatype=paramLoadrept['datatype'],
                        downloadonly=downloadonly, 
                        usePandas=False, usePyTplot=True)

    pytplot.del_data()
    varss_mageis = load_sat(trange=trange, satellite=paramLoadmageis['satellite'],
                        probe=paramLoadmageis['probe'], level=paramLoadmageis['level'], 
                        rel=paramLoadmageis['rel'], instrument=paramLoadmageis['instrument'],
                        datatype=paramLoadmageis['datatype'],
                        downloadonly=downloadonly, 
                        usePandas=False, usePyTplot=True)
    pytplot.del_data()
    varss_mageis = load_sat(trange=trange, satellite=paramLoadmageis['satellite'],
                        probe=paramLoadmageis['probe'], level='2', 
                        rel=paramLoadmageis['rel'], instrument=paramLoadmageis['instrument'],
                        datatype=paramLoadmageis['datatype'],
                        downloadonly=downloadonly, 
                        usePandas=False, usePyTplot=True)
    # quants_fedu_rept = pytplot.data_quants['FEDU']

    pytplot.del_data()
    varss_efw = load_sat(trange=trange, satellite=paramLoadEfw['satellite'],
                        probe=[probe], level=paramLoadEfw['level'], rel='rel03',
                        instrument=paramLoadEfw['instrument'], datatype=paramLoadEfw['datatype'],
                        varnames=paramLoadEfw['varnames'], downloadonly=downloadonly,
                        usePandas=False, usePyTplot=True)





#%%

trange= ['2012-08-01','2019-07-01']
config_file_sat = './pysatdata/resources/config_file.json'

varss_efw = downloadData(trange, config_file_sat)


#%%

# %%
