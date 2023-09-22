import pandas as pd
from pysatdata.utils.netcdf2tplot import *

#%%

def readData_goes(files, usePyTplot, usePandas, suffix, time):

    global time_ind, tvars
    for file in files:
        print(f'Reading... {file}')
        tvars = netcdf_to_tplot(file, suffix=suffix, merge=True, time=time)

    if usePyTplot:
        return tvars
    if usePandas:
        import pytz
        import datetime
        components = ["x", "y", "z"]
        out_dict = {}
        for vas in tvars:
            temp_var = pytplot.data_quants[vas].values
            if len(temp_var.shape) > 1:
                for i, c in enumerate(components):
                    out_dict[f"{vas}_{c}"] = list(temp_var[:, i])
            else:
                out_dict[f"{vas}"] = list(temp_var)

        for vas in tvars:
            if 'time' in pytplot.data_quants[vas].coords.keys():
                time_ind = pytplot.data_quants[vas].coords['time'].values
                break
        time = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time_ind]

        out_vars_df = pd.DataFrame(out_dict, index=time)

        return out_vars_df


