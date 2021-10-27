#%%
#
import gc
import pytz
from pysatdata.loaders.load import *
import datetime
from loguru import logger
from pysatdata.utils.flux2PhSD.flux2PhSD import *                       

#%%

def getData2PhSD(trange, dictParams, probeList, alphaRange, Kd, MUd):

    psdL = list()
    epochL = list()
    lstarL = list()
    for pp in probeList:
        pytplot.del_data()
        varss_reptEph = load_sat(trange=trange, satellite=dictParams['satellite'],
                            probe=[pp], level=dictParams['level'], rel=dictParams['rel']['ephem'],
                            instrument=dictParams['instrument']['ephem'],datatype=dictParams['datatype']['ephem'],
                            varnames=['Epoch', 'K', 'Alpha', 'Lstar', 'Lsimple', 'ApogeeTimes'],
                            config_file=config_file_sat, 
                            downloadonly=False,
                            usePandas=False, usePyTplot=True)

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
        # tempTimeList = list()
        # for i in apoT:
        #     tempTimeList.append(datetime.datetime.strptime(i, '%Y-%m-%dT%H:%M:%S.%fZ'))

        pytplot.del_data()
        varss_emfisis = load_sat(trange=trange, satellite=dictParams['satellite'],
                            probe=[pp], level=dictParams['level'], rel=dictParams['rel']['emfisis'],
                            instrument=dictParams['instrument']['emfisis'],datatype=dictParams['datatype']['emfisis'],
                            cadence='1sec', coord='gse',
                            varnames=['Magnitude'],
                            config_file=config_file_sat, downloadonly=False,
                            usePandas=False, usePyTplot=True)

        magVars = pytplot.data_quants['Magnitude']
        magB = magVars.values


        pytplot.del_data()
        varss_mageis = load_sat(trange=trange, satellite=dictParams['satellite'],
                            probe=[pp], level=dictParams['level'], rel=dictParams['rel']['mageis'],
                            instrument=dictParams['instrument']['mageis'],datatype=dictParams['datatype']['mageis'],
                            config_file=config_file_sat, 
                            downloadonly=False,
                            usePandas=False, usePyTplot=True)

        feduQuants = pytplot.data_quants['FEDU']
        feduMageis = feduQuants.values
        specbinMageis = feduQuants.coords['spec_bins'].values
        # feduAlpQuants = pytplot.data_quants['FEDU_Alpha']
        feduAlpha = feduQuants.coords['v1'].values
        timeMageis = feduQuants.coords['time'].values
        epochMageis = [datetime.datetime.fromtimestamp(i, pytz.timezone("UTC")) for i in timeMageis]
        feduMageis[feduMageis == -9999999848243207295109594873856.000] = np.nan
        feduMageis[feduMageis == -1e31] = np.nan

        
        fpd = flux2PhSD(alphaRange=alphaRange, Kd=Kd, MUd=MUd, epoch_eph=epoch, 
                        specbinMageis=specbinMageis, k_Ephem=Kvalues,  
                        lStar_eph=lstar, lShell_eph=lsimple, alpha_eph=alpha,
                        emfis_mag=magB, pAngle_mageis=feduAlpha, fedu_mageis=feduMageis)

        
        mud, kd = fpd.extractAndProcessData()

        
        fpd.calcPsd()

        psdL.append(fpd.psd)
        lstarL.append(fpd.Ls_forKd)
        epochL.append(fpd.Epo)

    return epochL, lstarL, psdL


#%%
dictParams = {"satellite": 'rbsp', 
                "level": '3', 
                "rel": {'ephem': "rel03",
                        'emfisis': "rel03",
                        'mageis': 'rel04'},
                "instrument": {'ephem':'ephemeris',
                               'emfisis': 'emfisis',
                               'mageis': 'mageis'}, 
                "datatype": {'ephem':'definitive',
                               'emfisis': 'magnetometer',
                               'mageis': 'sectors'}}

stringInstant = '2014-12-21 19:11:00'
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


config_file_sat = '/home/jose/python_projects/pySatData/pysatdata/resources/config_file.json'

alphaRange = [80,90]
Kd = None
alphaD = np.arange(alphaRange[0],alphaRange[1]) #desired value of pitch angle; change
# Kd=.1 #desired value of K; change
MUd=1500 #desired MU in MeV/G; change 

probeList=['a', 'b']

# %%
testes =  getData2PhSD(trange, dictParams, probeList, alphaRange, Kd, MUd)
# %%
