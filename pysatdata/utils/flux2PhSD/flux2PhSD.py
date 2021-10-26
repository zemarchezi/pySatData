#%%
import numpy as np
import glob, os
from decimal import Decimal
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import stats
from scipy.signal import find_peaks
from matplotlib import dates
from scipy.interpolate import interp1d
from scipy.optimize import leastsq
from scipy import stats
from time import mktime
import pandas as pd
import dateutil.parser
# from plt_functions import *
# rc('text', usetex=True)
#%%
class flux2PhSD():
    '''
    This class aims to convert electron flux measured by the Van Allen Probes 
    to phase space density values as a function of the three adiabatic invariants

    The ephemeris data used is the TS04 '*def_MagEphem_TS04D*' (availabe at: https://rbsp-ect.newmexicoconsortium.org/data_pub/rbspa/MagEphem/definitive/)
    EMFISIS data used is the rel04 '*magnetometer_1sec-gse_emfisis-l3*' (availabe at: https://cdaweb.gsfc.nasa.gov/pub/data/rbsp/)
    MagEis data used is the rel04 '*rel04_ect-mageis-L3*' (availabe at: https://cdaweb.gsfc.nasa.gov/pub/data/rbsp/)

    To run this you must have the data files for all days you want to analyse. The data must be stored in the right folder and the initial and final date must be given.

    Input:
       inidate -> Initial date (yyyy-mm-dd)
       enddate -> Final date (yyyy-mm-dd)
       alphaRange -> list with the initial and end of desired pitch angle range. i.e alphaRange = [50,60] 
       Kd -> desired value of K in G^0.5/G (If Kd = None, the value is calculated based on the alpha Range)
       MUd -> deside mu value (first adiabatic invariant, given in MeV/G)
       path_ephe -> path where ephemeris data are stored
       path_emfi -> path where EMFISIS instrument data are stored
       path_mageis -> path where Mageis instrument data are stored
       ylim -> list with y-axis limits on the plots.
       probe -> wich Van Allen Probe ('A' or 'B')

    
    References: 
    [1] Turner et al. 2012
    [2] Green & Kivelson, 2004.
    [3] Hartley & Denton, 2014  
    '''

    def __init__(self, alphaRange, Kd, MUd, path_ephe, path_emfi, path_mageis, ylim, probe):

        self.path_ephe = path_ephe
        self.path_emfi = path_emfi
        self.path_mageis = path_mageis
        self.alphaD = arange(alphaRange[0],alphaRange[1])
        self.Kd = Kd
        self.MUd = MUd
        self.ylim = ylim
        self.daysArray = pd.DatetimeIndex(arange(datetime64(inidate), datetime64(enddate)))
        self.probe = probe



    def extractAndProcessData(self):
        self.K_list = list()
        self.t_list = list()
        self.pAngle_list = list()
        self.lStar_list = list()
        self.b_emf_list = list()
        self.pAngle_mageis_list = list()
        self.fedu_mageis_list = list()
        self.epoch_mageis_list = list()
        self.lvalue_list = list()
        self.apogee_list = list()
        for dd in self.daysArray:
            #
            # data from the magnetic ephemeris
            #
            ff_eph = glob.glob(self.path_ephe + 'rbsp{:}_def_MagEphem_TS04D_{:}{:02d}{:02d}_*.h5'.format(probe.lower(), dd.year, dd.month, dd.day))
            eph_filename = ff_eph[0]

            data_eph = h5py.File(eph_filename, 'r') #vap 

            self.K_list.append(data_eph['K'][:])
            self.t_list.append(data_eph['IsoTime'][:])
            self.pAngle_list.append(data_eph['Alpha'][:])
            self.lStar_list.append(data_eph['Lstar'][:]) #Real lstar (magnetic field models older than TS04 might not resolve the adiabatic changes in L* as good as TS04) 
            self.lvalue_list.append(data_eph['Lsimple'][:])
            tempTimeList = list()
            for i in data_eph['ApogeeTimes'][:]:
                tempTimeList.append(pd.to_datetime(i.decode('utf-8')))
            self.apogee_list.append(tempTimeList)
            #
            # data from emfisis
            #
            ff_emfisis = glob.glob(self.path_emfi + 'rbsp-{:}_magnetometer_1sec-gse_emfisis-l3_{:}{:02d}{:02d}_*.cdf'.format(probe.lower(), dd.year, dd.month, dd.day))
            emfisis_filename = ff_emfisis[0]

            emfisis_data = pycdf.CDF(emfisis_filename)

            self.b_emf_list.append(emfisis_data['Magnitude'][:]) #real local B
            #
            # data from mageis
            #
            ff_mageis = glob.glob(self.path_mageis + 'rbsp{:}_rel04_ect-mageis-L3_{:}{:02d}{:02d}_*.cdf'.format(probe.lower(), dd.year, dd.month, dd.day))
            mageis_filename = ff_mageis[0]

            mageis_data = pycdf.CDF(mageis_filename)

            self.pAngle_mageis_list.append(mageis_data['FEDU_Alpha'][:]) #pa vector from mageis
            self.fedu_mageis_list.append(mageis_data['FEDU'][:])
            self.epoch_mageis_list.append(mageis_data['Epoch'][:])


        #%%
        apogeeTime = concatenate(tuple(self.apogee_list), axis=0)

        self.K = concatenate(tuple(self.K_list), axis=0)
        # create an relation of K and alpha
        kk = K
        kk[kk == -1.00000000e+31] = nan
        kk[kk < 0] = nan
        # get the mean K profile during the period
        self.asa = pd.DataFrame(K.T)
        self.mediaK = asa.mean(axis=1)

        self.numbPoints, self.numbPoints_angle = shape(K) #bins of time and PA
        # numbPoints = 2000

        time = concatenate(tuple(t_list), axis=0)
        self.ttime = []
        for i in time:
            self.ttime.append(dateutil.parser.parse(i))

        # ttime = stats.binned_statistic(arange(len(ttime)),ttime, 'mean', bins=numbPoints)[0]#

        self.pae = self.pAngle_list[0]

        self.Ls = concatenate(tuple(self.lStar_list), axis=0)
        self.Lval = concatenate(tuple(self.lvalue_list), axis=0)
        #####
        # data from van allen
        self.B = concatenate(tuple(self.b_emf_list), axis=0)
        self.B2=stats.binned_statistic(arange(len(self.B)),self.B, 'mean', bins=self.numbPoints)[0]#

        self.pAngle_mageis = concatenate(tuple(self.pAngle_mageis_list), axis=0)
        self.fedu_mageis = concatenate(tuple(self.fedu_mageis_list), axis=0)


        self.Epo = concatenate(tuple(self.epoch_mageis_list), axis=0)
        #%%


        if not self.Kd:
            self.Kd = kValueCalc(angleRange=self.alphaD,  K_Vector=self.mediaK, alpha_Vector=self.pae)

        if not self.MUd:
            self.MUd = muCalc()


        #%%
        a1,b1,c1=shape(self.fedu_mageis)#ntime,npa,nE
        self.Enm=(mageis_data['FEDU_Energy'][:])*1e-3 #mageis energy channels changed to MeV
        self.Enm=self.Enm[range(21)]#range of valid energy channels for mageis
        #
        self.fedu_m2=zeros((self.numbPoints,b1,c1))
        # Epo_m2 = zeros()
        for i in range(b1):
            for j in range(c1):
                temp=self.fedu_mageis[:,i,j]
                temp[temp==-1.00000000e+31] = nan
                temp[temp==-9999999848243207295109594873856.000]=nan
                self.fedu_m2[:,i,j]=stats.binned_statistic(arange(a1),self.fedu_mageis[:,i,j], 'mean', bins=self.numbPoints)[0]#

    def calcPsd(self):
        self.psd=zeros((self.numbPoints))
        self.Ls_K=zeros((self.numbPoints))
        self.Lval_K=zeros((self.numbPoints))#psd and L* final vectors
        self.E_mu_spec=zeros((self.numbPoints)) #energy radial coverage for the chosen MU
        for iit in range(self.numbPoints):
            #############################################################################
            #1ST STEP, getting alphaK and E_mu (energy respective to desired K and MU)
            K_temp=self.K[iit,:]
            K_temp[K_temp==-1.00000000e+31] = nan
            
            temp_nan0=isnan(K_temp)#removing nan values in K prior to the interpolation
            pos0=where(temp_nan0==False)[0]
            if len(pos0)<2 or self.Kd<nanmin(K_temp) or self.Kd>nanmax(K_temp):
                self.psd[iit]=nan
                self.Ls_K[iit]=nan
                self.Lval_K[iit]=nan
                self.E_mu_spec[iit]=nan
            #Kd not defined means that that population is in a region of B with opened field lines
            #For the same reason, Lstar is no longer defined
            else:    
                pa_st0=pos0[0]
                pa_end0=pos0[-1]+1
                K_temp2=K_temp[pa_st0:pa_end0]
                pae2=self.pae[pa_st0:pa_end0]
                
                fK = interp1d(K_temp2,pae2) #linear fit, VALIDATED!!!
                self.aK=fK(self.Kd) #Alpha for desired K at instant iit
                #disp(aK) ok, no nan values generated for aK during interpolation
                
                #################################################################
                #LSTAR step, Linear interpolation over L*(t,alpha) to get L* at desired K, that is, at calculated alphaK
                #With the assumption that K provided was calculated for local alpha
                #that is, along the satellite's orbit         
                Ls_temp=self.Ls[iit,:]        
                Ls_temp[Ls_temp==-1.00000000e+31] = nan
                
                temp_nan2=isnan(Ls_temp)#removing nan values of Ls vector before interpolation
                pos2=where(temp_nan2==False)[0]
                if len(pos2)<2:
                    self.Ls_K[iit]=nan
                    self.Lval_K[iit]=nan
                else:          
                    Ls_temp2=Ls_temp[pos2]
                    pae2=self.pae[pos2]
                    
                    fLs = interp1d(pae2,Ls_temp2) #linear fit on L*-->y in respect to PA-->x, 
                    if self.aK<min(pae2) or self.aK>max(pae2):#ok
                        self.Ls_K[iit]=nan
                    else:    
                        self.Ls_K[iit]=fLs(aK)
                
                m0=9.11e-31 # electron rest mass Kg
                c=3e8 #speed of light in m/s 
                #mc_sq=1e-6*(m0*c**2)/(1.6022e-19)=.511 MeV #electron energy rest in units of MeV
                
                p2=2*(.511/c**2)*(self.B2[iit]*1e-5)*self.MUd/(sin(radians(self.aK))**2) #relativistic momentum squared       
                #B is converted from nT to Gauss        
                delta=(2*.511)**2+(4*p2*c**2)
                #solving the second degree equation in respect to E to get the kinetic energy En_mu at desired MU and K         
                E_mu=(-2*.511+sqrt(delta))/2 #VALIDATED!!!
                self.E_mu_spec[iit]=E_mu        
                disp(E_mu)
                if E_mu<nanmin(Enm) or E_mu>nanmax(Enm):
                    self.psd[iit]=nan
                    #disp([E_mu,iit])
                else:    

                    #############################################################################
                    #2ND STEP, Linear fitting of the pitch angle distribution
                    #to get the flux at aK for different energy channels (all saved in jK)
                    jK=zeros((21))#12    
                    for en in range(21):#12
                        temp_nan=isnan(self.fedu_m2[iit,:,en])#finding nan values for flux in funtion of PA
                        pos=where(temp_nan==False)[0]
                        if len(pos)==0:
                            jK[en]=nan
                        else:    
                            pa_st=pos[0]
                            pa_end=pos[-1]+1
                            #disp([pa_st,pa_end])
                            x=self.pAngle_mageis[pa_st:pa_end]
                            y=self.fedu_m2[iit,pa_st:pa_end,en]
                            
                            if self.aK<x[0] or self.aK>x[-1]:
                                self.psd[iit]=nan
                            else:    
                                fPAD=interp1d(x,y)
                                jK[en]=fPAD(self.aK)
                        
                    ############################################################################
                    #3RD STEP, exponential fit between fluxes (jK) and the energy channels
                    #to get j(aK) at E_mu.
                    temp_nan1=isnan(jK)#removing nan values of jK vector before interpolation
                    pos1=where(temp_nan1==False)[0]
                    if len(pos1)==0:
                        psd[iit]=nan
                    else:    
                        jK2=jK[pos1]
                        Enm2=Enm[pos1]
                        
                        fE = interp1d(Enm2,log(jK2)) #linear fit on ln(jK)-->y in respect to E-->x, VALIDATED!!! 
                        if E_mu<min(Enm2) or E_mu>max(Enm2):#ok
                            self.psd[iit]=nan
                        else:                
                            lnj=fE(E_mu) #log(jK) for E_mu
                            jd=e**lnj#interpolated j at alphaK and E_mu
                            
                            #############################################################################
                            #4TH STEP, interpolated j is finally converted to psd 
                            psd[iit]=3.325*(1e-8)*(jd)/(E_mu*(E_mu+2*.511))
                            #equivalent to j at fixed K and MU divided by pÂ² (the relativistic momentum squared) 
                            #this psd is given in the GEM unit of (c/MeV/cm)Â³
                            #for mageis, jd does not need to be multiplied by 1e-3


#%%
# 
import gc

import pandas as pd
import pytz
from pysatdata.loaders.load import *
import datetime
from loguru import logger                            
# %%
stringInstant = '2017-05-28 12:00:00'
#
instDate = datetime.datetime.strptime(stringInstant.split(' ')[0], '%Y-%m-%d')


inidate = instDate - datetime.timedelta(days = 1)
enddate = instDate + datetime.timedelta(days = 2)
trange0 = [inidate.strftime('%Y-%m-%d'), enddate.strftime('%Y-%m-%d')]
trange= [(inidate - datetime.timedelta(days = 1)).strftime('%Y-%m-%d'),
         (enddate + datetime.timedelta(days = 3)).strftime('%Y-%m-%d')]

tt0 = datetime.datetime.strptime(stringInstant, '%Y-%m-%d %H:%M:%S') - datetime.timedelta(days = 1)
tt1 = datetime.datetime.strptime(stringInstant, '%Y-%m-%d %H:%M:%S') + datetime.timedelta(days = 2)
testMaxParametersIndex = [tt0, tt1]

# trange0 = ['2020-12-14', '2020-12-15']
# trange= ['2020-12-13', '2020-12-16']

config_file_sat = '/home/jose/python_projects/pySatData/pysatdata/resources/config_file.json'
#%%
varss_reptEph = load_sat(trange=trange, satellite='rbsp',
                     probe=['a'], level='3', rel='rel03',
                     instrument='ephemeris',datatype='definitive',
                     varnames=['Epoch', 'K', 'Alpha', 'Lstar', 'Lsimple', 'ApogeeTimes'],
                     config_file=config_file_sat, 
                     downloadonly=False,
                     usePandas=False, usePyTplot=True)
#%%
kK = pytplot.data_quants['K']
Kvalues = kK.values
Kvalues[Kvalues == -1.00000000e+31] = np.nan
Kvalues[Kvalues < 0] = np.nan
time = kK.coords['time'].values
alpha = kK.coords['spec_bins'].values[0]
epoch = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in time]
lstar = pytplot.data_quants['Lstar'].values
lsimple = pytplot.data_quants['Lsimple'].values
apoT = pytplot.data_quants['ApogeeTimes']['data'][0]
tempTimeList = list()
for i in apoT:
    tempTimeList.append(datetime.datetime.strptime(i, '%Y-%m-%dT%H:%M:%S.%fZ'))
#%%
pytplot.del_data()
varss_emfisis = load_sat(trange=trange, satellite='rbsp',
                     probe=['a'], level='3', rel='rel03',
                     instrument='emfisis',datatype='magnetometer',
                     cadence='1sec', coord='gse',
                     varnames=['Magnitude'],
                     config_file=config_file_sat, downloadonly=False,
                     usePandas=False, usePyTplot=True)
#%%
magVars = pytplot.data_quants['Magnitude']
magB = magVars.values
#%%

pytplot.del_data()
varss_mageis = load_sat(trange=trange, satellite='rbsp',
                     probe=['a'], level='3', rel='rel04',
                     instrument='mageis',datatype='sectors',
                     config_file=config_file_sat, 
                     downloadonly=False,
                     usePandas=False, usePyTplot=True)
#%%
feduQuants = pytplot.data_quants['FEDU']
feduMageis = feduQuants.values
# feduAlpQuants = pytplot.data_quants['FEDU_Alpha']
feduAlpha = feduQuants.coords['v1'].values
timeMageis = feduQuants.coords['time'].values
epoch = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in timeMageis]
feduMageis[feduMageis == -9999999848243207295109594873856.000] = np.nan
feduMageis[feduMageis == -1e31] = np.nan
#%%

mediaK = np.nanmean(Kvalues, axis=1)

numbPoints, numbPoints_angle = Kvalues.shape #bins of time and PA
# numbPoints = 2000
pae = alpha
Ls = lstar
Lval = lsimple
#%%

time = concatenate(tuple(t_list), axis=0)
ttime = []
for i in time:
    ttime.append(dateutil.parser.parse(i))

# ttime = stats.binned_statistic(arange(len(ttime)),ttime, 'mean', bins=numbPoints)[0]#



Ls = concatenate(tuple(lStar_list), axis=0)
Lval = concatenate(tuple(lvalue_list), axis=0)
#####
# data from van allen
B = concatenate(tuple(b_emf_list), axis=0)
B2=stats.binned_statistic(arange(len(B)),B, 'mean', bins=numbPoints)[0]#

pAngle_mageis = concatenate(tuple(pAngle_mageis_list), axis=0)
fedu_mageis = concatenate(tuple(fedu_mageis_list), axis=0)


Epo = concatenate(tuple(epoch_mageis_list), axis=0)


#%%


import cdflib
# %%
dd = cdflib.CDF('/home/jose/sat_data/rbsp/rbspa/l3/ect/mageis/sectors/rel03/2017/rbspa_rel03_ect-mageis-l3_20170529_v7.1.1.cdf')
# %%
