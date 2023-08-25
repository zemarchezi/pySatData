from pysatdata.utils.cdf_to_tplot import *
import warnings
import astropy
import pandas as pd

def readData_ace(files, usePyTplot, usePandas,
                  suffix, get_support_data,
                  varformat, varnames, notplot):

    global time_ind, tvars

    with warnings.catch_warnings():
        # for some reason, ACE CDFs throw ERFA warnings (likely while converting
        # times inside astropy); we're ignoring these here
        # see: https://github.com/astropy/astropy/issues/9603
        warnings.simplefilter('ignore', astropy.utils.exceptions.ErfaWarning)
        tvars = cdf_to_tplot(files, suffix=suffix, get_support_data=get_support_data,
                         varformat=varformat, varnames=varnames, notplot=notplot)

    

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
        out_vars_df = pd.DataFrame(out_dict)
        for vas in tvars:
            if 'time' in pytplot.data_quants[vas].coords.keys():
                time = pytplot.data_quants[vas].coords['time'].values
                # time = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time_ind]
                out_vars_df = pd.DataFrame(out_dict, index=time)
            else:
                out_vars_df = pd.DataFrame(out_dict)
        
    return out_vars_df
        

        
