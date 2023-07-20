import datetime
import numpy as np
import pandas as pd
from scipy import interpolate as interp
from pysatdata.utils.library_functions import fill_nan


def interpolateFluxRbsp(enSignal, lValues, timeArray, resolution_L=0.025):
    """
    enSignal: electron flux at one energy level (1D array)
    lValues: L-shell or l-star time series
    timeArray: datetime array 

    """

    L_inerp = fill_nan(lValues)

    tHour0 = timeArray[0].hour
    tHour1 = timeArray[-1].hour
    tMin0 = timeArray[0].minute
    tMin1 = timeArray[-1].minute
    tSec0 = timeArray[0].second
    tSec1 = timeArray[-1].second
    dad0 = int(timeArray[0].strftime('%j')) + ((tHour0 / 24) + (tMin0 / (60*24)) + (tSec0 / (60*60*24)))
    dad1 = int(timeArray[-1].strftime('%j')) + ((tHour1 / 24) + (tMin1 / (60*24)) + (tSec1 / (60*60*24)))

    time = np.linspace(dad0 , dad1, len(L_inerp))
    p = np.matrix.transpose(np.asmatrix([time, L_inerp]))
    z = enSignal

    xtime, xL = np.meshgrid(np.arange(min(time), max(time), 0.01), np.arange(min(L_inerp), max(L_inerp),resolution_L))

    # datetime x array
    xax = [datetime.datetime(timeArray[0].year, 1, 1, 0, 0) + datetime.timedelta(t - 1) for t in xtime[0,:]]

    # creade Y axis array (L values)
    yax = np.arange(min(L_inerp), max(L_inerp),resolution_L)

    # interpolate the data indo the grid
    flux = interp.griddata(p, z, (xtime, xL), method='linear')

    # Elimitade nan values
    flux[flux == 0.0] = 'NaN'
    maskflux = np.ma.masked_where(np.isnan(flux),flux, copy=True)

    return xax, yax, maskflux
