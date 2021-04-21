from pysatdata.loaders.load import *
import json

#%%

with open('/home/jose/python_projects/pySatData/pysatdata/utils/resources/config_file.json', 'r') as f:
    config_file_sat = json.load(f)


#%%
trange=['2013-04-12', '2013-04-19']
varss_rept = load_sat(trange=trange, satellite='rbsp',
                     probe=['a'], level='l3', rel='rel03',
                     instrument='rept',
                     config_file=config_file_sat, downloadonly=False,
                     usePandas=False, usePyTplot=True)
