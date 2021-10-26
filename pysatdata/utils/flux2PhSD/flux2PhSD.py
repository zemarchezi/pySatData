#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import stats
from scipy.interpolate import interp1d
from scipy import stats
import pandas as pd
import datetime
# from plt_functions import *
# rc('text', usetex=True)
#%%

def kValueCalc(angleRange, K_Vector, alpha_Vector):
    j_k = K_Vector
    k_alpha = alpha_Vector
    y = np.zeros(len(j_k))
    for k in range(len(j_k)):
        y[k]=np.log10(j_k[k])
    x = np.zeros(len(k_alpha))
    for k in range(len(k_alpha)):
        x[k]=k_alpha[k]
    f = interp1d(x, y, kind='linear',fill_value='extrapolate')


    calcK = list()
    for i in angleRange:
        calcK.append(10**f(i))

    return (np.mean(calcK))


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

    def __init__(self, alphaRange, Kd, MUd, epoch_eph, specbinMageis,
                 k_Ephem, lStar_eph, lShell_eph, alpha_eph,
                 emfis_mag, pAngle_mageis, fedu_mageis):

        self.Kvalues = k_Ephem
        self.Ls = lStar_eph
        self.Lval = lShell_eph
        self.Epo = epoch_eph
        self.pae = alpha_eph
        self.B = emfis_mag
        self.pAngle_mageis = pAngle_mageis
        self.fedu_mageis = fedu_mageis
        self.specbinMageis = specbinMageis
        self.alphaD = np.arange(alphaRange[0],alphaRange[1])
        self.Kd = Kd
        self.MUd = MUd
        



    def extractAndProcessData(self):
        

        self.mediaK = np.nanmean(self.Kvalues, axis=0)

        self.numbPoints, self.numbPoints_angle = self.Kvalues.shape #bins of time and PA
        # numbPoints = 2000

        # ttime = stats.binned_statistic(arange(len(ttime)),ttime, 'mean', bins=numbPoints)[0]#

        # self.lvalue_list
        #####
        # data from van allen

        self.B2=stats.binned_statistic(np.arange(len(self.B)),self.B, 'mean', bins=self.numbPoints)[0]#



        if not self.Kd:
            self.Kd = kValueCalc(angleRange=self.alphaD,  K_Vector=self.mediaK, alpha_Vector=self.pae)

        if not self.MUd:
            self.MUd = muCalc()



        a1,b1,c1=self.fedu_mageis.shape #ntime,npa,nE
        self.Enm=(self.specbinMageis)*1e-3 #mageis energy channels changed to MeV
        self.Enm=self.Enm[range(21)]#range of valid energy channels for mageis
        #
        self.fedu_m2=np.zeros((self.numbPoints,b1,c1))
        # Epo_m2 = zeros()
        for i in range(b1):
            for j in range(c1):
                temp=self.fedu_mageis[:,i,j]
                temp[temp==-1.00000000e+31] = 0
                temp[temp==-9999999848243207295109594873856.000]=0
                self.fedu_m2[:,i,j]=stats.binned_statistic(np.arange(a1),self.fedu_mageis[:,i,j], 'mean', bins=self.numbPoints)[0]#
        
        return (self.MUd, self.Kd)
    def calcPsd(self):
        self.psd=np.zeros((self.numbPoints))
        self.Ls_K=np.zeros((self.numbPoints))
        self.Lval_K=np.zeros((self.numbPoints))#psd and L* final vectors
        self.E_mu_spec=np.zeros((self.numbPoints)) #energy radial coverage for the chosen MU
        for iit in range(self.numbPoints):
            #############################################################################
            #1ST STEP, getting alphaK and E_mu (energy respective to desired K and MU)
            K_temp=self.Kvalues[iit,:]
            K_temp[K_temp==-1.00000000e+31] = np.nan
            
            temp_nan0=np.isnan(K_temp)#removing nan values in K prior to the interpolation
            pos0=np.where(temp_nan0==False)[0]
            if len(pos0)<2 or self.Kd<np.nanmin(K_temp) or self.Kd>np.nanmax(K_temp):
                self.psd[iit]=np.nan
                self.Ls_K[iit]=np.nan
                self.Lval_K[iit]=np.nan
                self.E_mu_spec[iit]=np.nan
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
                # Ls_temp=self.Ls#[iit,:]   
                Ls_temp=self.Ls[iit,:]        
                Ls_temp[Ls_temp==-1.00000000e+31] = np.nan
                
                temp_nan2=np.isnan(Ls_temp)#removing nan values of Ls vector before interpolation
                pos2=np.where(temp_nan2==False)[0]
                if len(pos2)<2:
                    self.Ls_K[iit]=np.nan
                    self.Lval_K[iit]=np.nan
                else:          
                    Ls_temp2=Ls_temp[pos2]
                    pae2=self.pae[pos2]
                    
                    fLs = interp1d(pae2,Ls_temp2) #linear fit on L*-->y in respect to PA-->x, 
                    if self.aK<min(pae2) or self.aK>max(pae2):#ok
                        self.Ls_K[iit]=np.nan
                    else:    
                        self.Ls_K[iit]=fLs(self.aK)
                
                m0=9.11e-31 # electron rest mass Kg
                c=3e8 #speed of light in m/s 
                #mc_sq=1e-6*(m0*c**2)/(1.6022e-19)=.511 MeV #electron energy rest in units of MeV
                
                p2=2*(.511/c**2)*(self.B2[iit]*1e-5)*self.MUd/(np.sin(np.radians(self.aK))**2) #relativistic momentum squared       
                #B is converted from nT to Gauss        
                delta=(2*.511)**2+(4*p2*c**2)
                #solving the second degree equation in respect to E to get the kinetic energy En_mu at desired MU and K         
                E_mu=(-2*.511+np.sqrt(delta))/2 #VALIDATED!!!
                self.E_mu_spec[iit]=E_mu        
                np.disp(E_mu)
                if E_mu<np.nanmin(self.Enm) or E_mu>np.nanmax(self.Enm):
                    self.psd[iit]=np.nan
                    #disp([E_mu,iit])
                else:    

                    #############################################################################
                    #2ND STEP, Linear fitting of the pitch angle distribution
                    #to get the flux at aK for different energy channels (all saved in jK)
                    jK=np.zeros((21))#12    
                    for en in range(21):#12
                        temp_nan=np.isnan(self.fedu_m2[iit,:,en])#finding nan values for flux in funtion of PA
                        pos=np.where(temp_nan==False)[0]
                        if len(pos)==0:
                            jK[en]=np.nan
                        else:    
                            pa_st=pos[0]
                            pa_end=pos[-1]+1
                            #disp([pa_st,pa_end])
                            x=self.pAngle_mageis[pa_st:pa_end]
                            y=self.fedu_m2[iit,pa_st:pa_end,en]
                            
                            if self.aK<x[0] or self.aK>x[-1]:
                                self.psd[iit]=np.nan
                            else:    
                                fPAD=interp1d(x,y)
                                jK[en]=fPAD(self.aK)
                        
                    ############################################################################
                    #3RD STEP, exponential fit between fluxes (jK) and the energy channels
                    #to get j(aK) at E_mu.
                    temp_nan1=np.isnan(jK)#removing nan values of jK vector before interpolation
                    pos1=np.where(temp_nan1==False)[0]
                    if len(pos1)==0:
                        self.psd[iit]=np.nan
                    else:    
                        jK2=jK[pos1]
                        Enm2=self.Enm[pos1]
                        
                        fE = interp1d(Enm2,np.log(jK2)) #linear fit on ln(jK)-->y in respect to E-->x, VALIDATED!!! 
                        if E_mu<min(Enm2) or E_mu>max(Enm2):#ok
                            self.psd[iit]=np.nan
                        else:                
                            lnj=fE(E_mu) #log(jK) for E_mu
                            jd=np.e**lnj#interpolated j at alphaK and E_mu
                            
                            #############################################################################
                            #4TH STEP, interpolated j is finally converted to psd 
                            self.psd[iit]=3.325*(1e-8)*(jd)/(E_mu*(E_mu+2*.511))
                            #equivalent to j at fixed K and MU divided by pÂ² (the relativistic momentum squared) 
                            #this psd is given in the GEM unit of (c/MeV/cm)Â³
                            #for mageis, jd does not need to be multiplied by 1e-3

    
    def separateOrbits(self, ylim):
        self.ylim = ylim
        self.dfPsd = pd.DataFrame(self.psd,index=self.ttime)
        # disp(self.dfPsd)
        self.dfPsd['lStar'] = self.Ls_K

        dfP = self.dfPsd.copy()
        dfP = dfP.interpolate()
        dfP[dfP['lStar'] <=1.5] = np.nan
        # dfP[dfP['lStar'] >=5.2] = np.nan
        # dfP[0] = self.dfPsd[0]

        orbits = np.split(dfP, np.where(np.isnan(dfP['lStar'].values))[0])
        # removing NaN entries
        orbits = [ev[~np.isnan(ev.values)] for ev in orbits if not isinstance(ev, np.ndarray)]
        # removing empty DataFrames
        orbits = [ev for ev in orbits if not ev.empty]

        self.phDenInB = list()
        self.phDenOuB = list()
        for nn, i in enumerate(orbits):
            maximIndex = i['lStar'].idxmax()
            slice_out = i.loc[i.index[0]:maximIndex]
            # slice_out[slice_out['lStar'] >=5.6] = np.nan
            slice_in = i.loc[maximIndex:i.index[-1]]
            # slice_in[slice_in['lStar'] >=5.6] = np.nan
            # self.phDenOuB.append(slice_out)
            # self.phDenInB.append(slice_in)
            self.phDenOuB.append(slice_out.resample('10T').mean())
            self.phDenInB.append(slice_in.resample('10T').mean())

        return (self.phDenInB, self.phDenOuB)


