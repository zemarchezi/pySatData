from scipy.signal import butter, lfilter, filtfilt
from scipy import interpolate as interp
import numpy as np
from numpy import pi, cos, sin, arctan2, sqrt, dot
import datetime
import pandas as pd
import aacgmv2
import igrf12
import pyIGRF


def normD(a):
    norm = 0
    for i in range(3):
        norm += a[i] * a[i]
    return np.sqrt(norm)


def crossD(a, b):
    cross = [0] * 3
    cross[0] = a[1] * b[2] + a[2] * b[1]
    cross[1] = a[2] * b[0] + a[0] * b[2]
    cross[2] = a[0] * b[1] + a[1] * b[2]
    return cross


def replace_at_index1(tup, ix, val):
    lst = list(tup)
    for i in range(0, len(ix)):
        lst[ix[i]] = val[i]
    return tuple(lst)


# Butterword filter coefficients
def butter_bandpass(lowcut, highcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


# Band-pass butterword filter
def butter_bandpass_filter(data, lowcut, highcut, fs, order):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


# fill the gaps with nans
def fill_nan(A):
    '''
     interpolate to fill nan values
     '''
    if np.isnan(A[0]):
        A[0] = 0.5
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = interp.interp1d(inds[good], A[good], bounds_error=False)
    B = np.where(np.isfinite(A), A, f(inds))
    return B


# convert the coordinate system from gse to field aligned system
def rotate_field_fac(x, y, z, bx, by, bz, ex, ey, ez):
    '''
    rotate the fields into Field Alignet Coordinate System

    data: pandas dataframe with the columns: 'x', 'y', 'ex', 'ey', 'ez', 'bx', 'by', 'bz'
    '''
    x = x.values
    y = y.values
    z = z.values
    v1x = np.transpose(ex)
    v1y = ey
    v1z = ez
    bx = bx.values
    by = by.values
    bz = bz.values

    # v1p = v1a = v1r = bp = ba = br = r =  b_fac = b_orig = np.zeros((len(x)))
    tempB = np.zeros((len(x), 3))
    tempE = np.zeros((len(x), 3))
    for i in range(0, len(x)):
        Jac = np.zeros((3, 3))
        r = [x[i], y[i], z[i]] / np.sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i])
        Jac[0, :] = [bx[i], by[i], bz[i]] / np.sqrt(bx[i] * bx[i] + by[i] * by[i] + bz[i] * bz[i])
        Jac[1, :] = crossD(Jac[0, :], r) / normD(crossD(Jac[0, :], r))
        Jac[2, :] = crossD(Jac[1, :], Jac[0, :]) / normD(crossD(Jac[1, :], Jac[0, :]))
        # apply rotation for B
        tempB[i, :] = np.dot(Jac, ([bx[i], by[i], bz[i]]))
        # Apply the rotation for vector V1
        tempE[i, :] = np.dot(Jac, ([v1x[i], v1y[i], v1z[i]]))
    #        # testing whether the rotation is correct
    #        tempFields['b_fac'][i] = np.linalg.norm([tempFields['bp'][i], tempFields['ba'][i], tempFields['br'][i]])
    #        tempFields['b_orig'][i] = np.linalg.norm([bx[i], by[i], bz[i]])

    return (pd.DataFrame(np.transpose([tempB[:, 0], tempB[:, 1], tempB[:, 2], tempE[:, 0], tempE[:, 1], tempE[:, 2]]),
                         columns=['bp', 'ba', 'br', 'v1p', 'v1a', 'v1r']))


def l_dipole(cgm_lat):
    return 1. / (np.cos(np.deg2rad(cgm_lat)) ** 2.)


def geo2mag(incoord):
    """geographic coordinate to magnetic coordinate:

        incoord is numpy array of shape (2,*)
        array([[glat0,glat1,glat2,...],
            [glon0,glon1,glon2,...])
        where glat, glon are geographic latitude and longitude
        (or if you have only one point it is [[glat,glon]])

        returns
        array([mlat0,mlat1,...],
            [mlon0,mlon1,...]])
        """

    # SOME 'constants'...
    lon = 360 - 72.6  # or 71.41W
    lat = 80.4
    r = 1.0

    # convert first to radians
    lon, lat = [x * pi / 180 for x in (lon, lat)]

    glat = incoord[0] * pi / 180.0
    glon = incoord[1] * pi / 180.0
    galt = glat * 0. + r

    coord = np.vstack([glat, glon, galt])

    # convert to rectangular coordinates
    x = coord[2] * cos(coord[0]) * cos(coord[1])
    y = coord[2] * cos(coord[0]) * sin(coord[1])
    z = coord[2] * sin(coord[0])
    xyz = np.vstack((x, y, z))

    # computer 1st rotation matrix:
    geo2maglon = np.zeros((3, 3), dtype='float64')
    geo2maglon[0, 0] = cos(lon)
    geo2maglon[0, 1] = sin(lon)
    geo2maglon[1, 0] = -sin(lon)
    geo2maglon[1, 1] = cos(lon)
    geo2maglon[2, 2] = 1.
    out = dot(geo2maglon, xyz)

    tomaglat = np.zeros((3, 3), dtype='float64')
    tomaglat[0, 0] = cos(.5 * pi - lat)
    tomaglat[0, 2] = -sin(.5 * pi - lat)
    tomaglat[2, 0] = sin(.5 * pi - lat)
    tomaglat[2, 2] = cos(.5 * pi - lat)
    tomaglat[1, 1] = 1.
    out = dot(tomaglat, out)

    mlat = arctan2(out[2],
                   sqrt(out[0] * out[0] + out[1] * out[1]))
    mlat = mlat * 180 / pi
    mlon = arctan2(out[1], out[0])
    mlon = mlon * 180 / pi

    # outcoord = np.vstack((mlat, mlon))
    return [mlat, (360.0 + mlon)]

#%%
def convert_coords(date: str = '20210405',
                   lat_long: list = [],
                   altitude_km: float = 100.):
    dt = datetime.datetime.strptime(date, '%Y%m%d')

    if dt.year <=2020:
        mag = igrf12.igrf(dt, glat=float(lat_long[0]), glon=float(lat_long[1]), alt_km=altitude_km)
        decl = mag['decl'].values[0]
    else:
        mag = pyIGRF.igrf_value(lat=float(lat_long[0]), lon=float(lat_long[1]), alt=altitude_km, year=dt.year)
        mag = {'decl': mag[0],
               'incl': mag[1],
               'horiz': mag[2],
               'north': mag[3],
               'east': mag[4],
               'down': mag[5],
               'total': mag[6]
        }
        decl = mag['decl']


    # igrf13: d, i, h, x, y, z, f

    cgm_lat, cgm_lon, cgm_r = aacgmv2.convert_latlon(lat_long[0],
                                                     lat_long[1],
                                                     altitude_km, dt, method_code='G2A')

    if cgm_lon < 0:
        cgm_lon = 360 + cgm_lon

    mlt_ut = aacgmv2.convert_mlt(cgm_lon, datetime.datetime(int(dt.year), 1, 1, 0, 0, 0), m2a=False)
    mlt = 24 - mlt_ut

    l_dip = 1. / (np.cos(np.deg2rad(cgm_lat)) ** 2.)

    out_dict = {
        'cgm_latitude': cgm_lat,
        'cgm_longitude': cgm_lon,
        'declination': decl,
        'lshell': l_dip,
        'mlt_midnight': mlt[0],
        'mlt_ut': mlt_ut[0]
    }

    return out_dict
