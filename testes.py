# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

#%%
from pysatdata.loaders.load import *
import json

with open('./pysatdata/utils/resources/config_file.json', 'r') as f:
    config_file = json.load(f)

#%%

trange=['2021-04-05', '2021-04-12']
varss = load(trange=trange, satellite='goes',
             probe=[16], level=2,
             instrument='magn', datatype='hires',
             config_file=config_file, downloadonly=False)



#%%

