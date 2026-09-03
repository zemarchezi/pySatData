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

    L_interp = np.asarray(fill_nan(lValues), dtype=float)
    z = np.asarray(enSignal, dtype=float)

    tHour0, tMin0, tSec0 = timeArray[0].hour, timeArray[0].minute, timeArray[0].second
    tHour1, tMin1, tSec1 = timeArray[-1].hour, timeArray[-1].minute, timeArray[-1].second
    dad0 = int(timeArray[0].strftime('%j')) + (tHour0/24 + tMin0/(60*24) + tSec0/(60*60*24))
    dad1 = int(timeArray[-1].strftime('%j')) + (tHour1/24 + tMin1/(60*24) + tSec1/(60*60*24))

    time = np.linspace(dad0, dad1, len(L_interp))

    # --- ponto-chave da correção ---
    # remove qualquer amostra onde time, L ou o próprio flux não sejam finitos
    valid = np.isfinite(time) & np.isfinite(L_interp) & np.isfinite(z)
    n_dropped = (~valid).sum()
    if n_dropped:
        print(f"Descartando {n_dropped} pontos inválidos (NaN/inf) antes da interpolação")

    if valid.sum() < 4:
        raise ValueError("Pontos válidos insuficientes para interpolar (restaram "
                          f"{valid.sum()} de {len(valid)}).")

    time_v = time[valid]
    L_v = L_interp[valid]
    z_v = z[valid]

    p = np.column_stack((time_v, L_v))

    xtime, xL = np.meshgrid(
        np.arange(time_v.min(), time_v.max(), 0.01),
        np.arange(L_v.min(), L_v.max(), resolution_L)
    )

    # eixo X em datetime
    xax = [datetime.datetime(timeArray[0].year, 1, 1, 0, 0) + datetime.timedelta(t - 1)
           for t in xtime[0, :]]

    # eixo Y (valores de L)
    yax = np.arange(L_v.min(), L_v.max(), resolution_L)

    flux = interp.griddata(p, z_v, (xtime, xL), method='linear')

    flux[flux == 0.0] = np.nan
    maskflux = np.ma.masked_where(np.isnan(flux), flux, copy=True)

    return xax, yax, maskflux