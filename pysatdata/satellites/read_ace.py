from pytplot import cdf_to_tplot

def readData_ace(files, usePyTplot, usePandas,
                  suffix, get_support_data,
                  varformat, varnames, notplot):

    global time_ind, tvars

    tvars = cdf_to_tplot(files, suffix=suffix, get_support_data=get_support_data,
                         varformat=varformat, varnames=varnames, notplot=notplot)

    if usePyTplot:
        return tvars
    if usePandas:
        pass
