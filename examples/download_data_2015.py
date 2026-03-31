#%%
from pysatdata.loaders.load import *
from pysatdata.utils.interpolate_flux_rbsp import * 
from pysatdata.utils.plotFunc.plot_funtions import *
from pysatdata.utils.library_functions import *
#%%
paramLoadrept = {"satellite": 'rbsp', "probe": ['b'], "level": '3', "rel": "rel03",
                    "instrument": 'rept', "datatype": "sectors"}
trange= ['2015-01-01','2016-01-01']
downloadonly = True

varss_rept = load_sat(trange=trange, satellite=paramLoadrept['satellite'],
                        probe=paramLoadrept['probe'], level=paramLoadrept['level'], 
                        rel=paramLoadrept['rel'], instrument=paramLoadrept['instrument'],
                        datatype=paramLoadrept['datatype'],testRemotePath=False,
                        downloadonly=downloadonly, 
                        searchFilesFirst=True, usePandas=False, usePyTplot=True)

