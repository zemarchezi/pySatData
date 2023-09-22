from scipy.signal import butter, lfilter, filtfilt
from scipy import interpolate as interp
import numpy as np
from numpy import pi, cos, sin, arctan2, sqrt, dot
import matplotlib.dates as mdates
import datetime
import math
import pandas as pd
from loguru import logger as logging
import urllib3
from urllib3 import PoolManager
import glob
#import numba
#from numba import jit

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

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

# fill the gaps with nans
def fill_nan(A):
    '''
     interpolate to fill nan values
     '''
    if np.isnan(A[0]):
        A[0] = np.nanmin(A)
    if np.isnan(A[-1]):
        A[-1] = np.nanmin(A)
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = interp.interp1d(inds[good], A[good], bounds_error=False)
    B = np.where(np.isfinite(A), A, f(inds))
    return B

def format_func(value, tick_number):
    hora = pd.to_datetime(mdates.num2date(value)).replace(tzinfo=None)

    return ('{:02d}:{:02d} UTC \n {:04d}/{:02d}/{:02d}'.format(hora.hour, hora.minute, hora.year, hora.month, hora.day))

#@numba.jit(nopython=True, nogil=True)
def calcExEFW(efield, bfield, filter=True):
    ey = efield[:,1]
    ez = efield[:,2]
    bx = bfield[:,0]
    by = bfield[:,1]
    bz = bfield[:,2]

    ex = np.zeros((len(bx)))
    for i in range(len(bx)):
        if filter:
            angle = math.degrees(math.atan(bx[i]/(math.sqrt(by[i]**2 + bz[i]**2))))
            if angle < 6.0:
                ex[i] = np.nan
            else:
                ex[i] = -((ey[i]*by[i]) + (ez[i]*bz[i]))/(bx[i])
        else:
            ex[i] = -((ey[i]*by[i]) + (ez[i]*bz[i]))/(bx[i])
    return ex
#%%
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
# convert the coordinate system from gse to field aligned system
def rotate_field_fac(x, y, z, bx, by, bz, ex, ey, ez):
    '''
    rotate the fields into Field Alignet Coordinate System

    data: pandas dataframe with the columns: 'x', 'y', 'ex', 'ey', 'ez', 'bx', 'by', 'bz'
    '''
    # v1p = v1a = v1r = bp = ba = br = r =  b_fac = b_orig = np.zeros((len(x)))
    bxs = smooth(bx,11)
    bys = smooth(by,11)
    bzs = smooth(bz,11)
    exs = smooth(ex,11)
    eys = smooth(ey,11)
    ezs = smooth(ez,11)
    tempB = np.zeros((len(x), 3))
    tempE = np.zeros((len(x), 3))
    for i in range(0, len(x)):
        JacB = np.zeros((3, 3))
        JacE = np.zeros((3, 3))
        r = [x[i], y[i], z[i]] / np.sqrt(x[i] **2 + y[i] **2 + z[i] **2)
        JacB[0, :] = [bxs[i], bys[i], bzs[i]] / np.sqrt(bxs[i] **2 + bys[i] ** 2 + bzs[i] **2 )
        JacB[1, :] = crossD(JacB[0, :], r) / normD(crossD(JacB[0, :], r))
        JacB[2, :] = crossD(JacB[1, :], JacB[0, :]) / normD(crossD(JacB[1, :], JacB[0, :]))
        JacE[0, :] = [exs[i], eys[i], ezs[i]] / np.sqrt(exs[i] **2 + eys[i] ** 2 + ezs[i] **2 )
        JacE[1, :] = crossD(JacE[0, :], r) / normD(crossD(JacE[0, :], r))
        JacE[2, :] = crossD(JacE[1, :], JacE[0, :]) / normD(crossD(JacE[1, :], JacE[0, :]))
        # apply rotation for B
        tempB[i, :] = np.dot(JacB, ([bx[i], by[i], bz[i]]))
        # Apply the rotation for E
        tempE[i, :] = np.dot(JacE, ([ex[i], ey[i], ez[i]]))

    temp_data = [tempB[:, 0], tempB[:, 1], tempB[:, 2], tempE[:, 0], tempE[:, 1], tempE[:, 2], x, y, z]
    return (pd.DataFrame(np.transpose(temp_data),
                         columns=['bp', 'ba', 'br', 'ep', 'ea', 'er','x', 'y', 'z']))


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


def cutFlux_lshell(enSignal,lValue, EnChanel, lArray, timeArray):

    l = float(lValue)
    cutF = list()
    cut_date = list()
    for i, ll in enumerate(lArray):
        if ll > l-0.01 and ll < l+0.01:
            cutF.append(enSignal[i, EnChanel])
            cut_date.append(timeArray[i])

    return cut_date, cutF

def calc_fce(bMagnitude):

    q_e = 1.60217653e-19
    m_e = 9.1094e-31 #electron's rest mass (kg)
    wce = q_e*(bMagnitude*1e-9)/m_e #Electron cyclotron frequency (rad/s)
    fce_calc = wce/(2*np.pi) #Electron cyclotron frequency (Hz)
    return fce_calc


def cutFlux_lshell2(enSignal, lvalue):
    cutF = enSignal.copy()
    l = float(lvalue)
    mask = (cutF['L'] < l-0.01)
    cutF[mask] = np.nan
    mask = (cutF['L'] > l + 0.01)
    cutF[mask] = np.nan

    cutF[cutF < 1e-10] = np.nan



    return cutF.interpolate('linear')

def testFiles(local_path, remote_names):
    non_exiting_names = []
    existing_names = []
    for rn in remote_names:
        temp_find = glob.glob(f"{local_path}/{rn}")
        if len(temp_find)==0:
            non_exiting_names.append(rn)
        else:
            existing_names.append(temp_find[0])
            logging.info(f"File is current: {temp_find[0]}")


    return non_exiting_names, existing_names


def testRemoteDir(config_file, satellite, prb, instrument, level, datatype):
    pool = PoolManager()
    logging.info("Testing Connection")
    try:
        remote_path = config_file[satellite]['remote_data_dir']
        responseSubpath = pool.request("GET", remote_path, preload_content=False,
                                    timeout=1)
        responseSubpath.close()
        logging.warning(f"Using {config_file[satellite]['remote_data_dir']}...")
        changeRemoteDir = False
    except (Exception) as e:
        # logger.error(e.args)
        logging.error(f"Directory {config_file[satellite]['remote_data_dir']} is not available...")
        remote_path = config_file[satellite]['remote_subpath'][str(prb)][instrument]["secondRemoteDir"]
        logging.info(f"testing {remote_path}...")

        responseSubpath = pool.request("GET", remote_path, preload_content=False,
                                    timeout=1)
        responseSubpath.close()
        logging.warning(f"Using {remote_path}...")
        changeRemoteDir = True
    except (Exception) as e:
        logging.error(e)
        logging.warning(f"There are no repository available")
        changeRemoteDir = False
        raise

    if changeRemoteDir:
        subpathKey = "altern_subpath"
        filenameKey = "altern_filename"
        if instrument in ['rept', 'mageis'] and level == '3':
            datatype = 'pitchangle'
    else:
        subpathKey = 'subpath'
        filenameKey = 'filename'
    return remote_path, subpathKey, filenameKey, datatype