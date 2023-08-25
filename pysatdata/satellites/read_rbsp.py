from pysatdata.utils.cdf_to_tplot import *

def readData_rbsp(files, usePyTplot, usePandas,
                  suffix, get_support_data,
                  varformat, varnames, notplot):

    global time_ind, tvars

    tvars = cdf_to_tplot(files, suffix=suffix, get_support_data=get_support_data,
                         varformat=varformat, varnames=varnames, notplot=notplot)

    if usePyTplot:
        return tvars
    if usePandas:
        import pandas as pd
        import pytz
        import datetime
        components = ["x", "y", "z"]
        out_dict = {}
        for vas in tvars:
            # print(vas)
            temp_var = pytplot.data_quants[vas].values
            if len(temp_var.shape) > 1:
                for i, c in enumerate(components):
                    out_dict[f"{vas}_{c}"] = list(temp_var[:, i])
            else:
                out_dict[f"{vas}"] = list(temp_var)

        for vas in tvars:
            if 'time' in pytplot.data_quants[vas].coords.keys():
                time = pytplot.data_quants[vas].coords['time'].values
                break
        # time = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time_ind]

        out_vars_df = pd.DataFrame(out_dict, index=time)

        return out_vars_df
