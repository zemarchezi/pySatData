#%%
from pysatdata.utils.cdf_to_tplot import *
import pandas as pd
import numpy as np
import datetime
import gzip

#%%
def readData_ace_merged(files, usePyTplot, usePandas,
                  suffix, get_support_data,
                  varformat, varnames, notplot):

    global time_ind, tvars
    columns = ["Year", "Decimal Day", "Hour", "Minute",
           "X_GSE", "Y_GSE", "Z_GSE", "B_Magnitude", 
           "Bx_GSE", "By_GSE", "Bz_GSE", "No_of_Vectors",
           "Quality_flag", "Proton_density", "Temperature",
           "HeH_Ratio", "Flow_speed", "Vx", "Vy", "Vz"]
    al_data = []
    for ff in files:
        da = []
        with gzip.open(ff,'rt') as f:
            for line in f:
                da.append(np.asarray(list(filter(None, line.split("\n")[0].split(" ")))).astype(float))

        da = np.asarray(da)
        timeDt = [datetime.datetime.strptime(f'{int(da[i][0])}{int(da[i][1])} {int(da[i][2])}:{int(da[i][3])}', '%Y%j %H:%M') for i in range(len(da))]

        mergedData = pd.DataFrame(da, index=timeDt, columns=columns)
        
        al_data.append(mergedData)
    
    mergedAce = pd.concat(al_data)

    if usePyTplot:
        print("Not implemented yet, returning pandas Dataframe")
        return mergedAce
    if usePandas:
        return mergedAce
        
        

        

# %%
