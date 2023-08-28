# Basic documentation for pySatData

Most functions' structure and some loaders are from the pySPEDAS (https://github.com/spedas/pyspedas) repository. I just used and organized it in a way that I felt was more accessible for my purposes :). There is an option to store the variables in a pandas DataFrame format or pytplot format.

## Requirements

Python 3.8+
### Installation: 
```
conda create -n yourenvname python=x.x anaconda
conda activate yourenvname

pip install pysatdata
```

***
## config_file.json

It can be done in this file if some change in downloading directories is needed.

```pysatdata/resources/config_file.json```

This file sets the HTTP directory for downloading data and the local directory to save the downloaded data.

The local directory files are organized as http directory, i.e.:
* RBSP: HOME/sat_data/rbsp/rbspa/l2/ect/rept/sectors/rel03/YYYY/filename.cdf

This files handles with the different http subpaths for the different levels and intruments in each probe.

***

## Examples: Loading data.

### RBSP REPT data
```python
#Import the loading functions.
from pysatdata.loaders.load import *

# Define the time range for downloading data.
trange=['2021-05-26', '2021-05-30']

# Loading Van Allen probes REPT data.
varss_rept = load_sat(trange=trange, satellite="rbsp",
                        probe=['a', 'b'], level="3", 
                        rel="rel03", instrument="rept",
                        datatype="sectors",
                        downloadonly=False, 
                        searchFilesFirst=True,
                        usePandas=False,
                        usePyTplot=True)
```
### RBSP EMFISIS data.

```python
varss_emfisis = load_sat(trange=trange, satellite="rbsp",
                            probe=['a','b'], rel="rel03", level="3",
                            instrument="emfisis", datatype="magnetometer",
                            cadence="1sec", coord="gse",
                            varnames=[], downloadonly=False,
                            usePandas=True, usePyTplot=False)
```

### RBSP MAGEIS data.

```python
varss_mageis = load_sat(trange=trange, satellite="rbsp",
                        probe=['a','b'], level="3", 
                        rel="rel03", instrument="mageis",
                        datatype="sectors",
                        downloadonly=False, 
                        usePandas=False, usePyTplot=True)
```
### RBSP EFW data.

```python
varss_efw = load_sat(trange=trange, satellite="rbps",
                        probe=['a','b'], level="2", rel='rel03',
                        instrument="efw", datatype="esvy_despun",
                        varnames=['efield_mgse', 'lshell'], downloadonly=False,
                        usePandas=False, usePyTplot=True)
```

### OMNI Solar wind data

```python
varss_omni = load_sat(trange=trange, satellite="omni",
                         probe="omni"
                         instrument="omni_cdaweb",datatype="hro_1min",
                         downloadonly=False,
                         usePandas=False, usePyTplot=True)
```

See ```plot_interpFlux_RBSP.py``` for an example of plotting the interpolated electron flux.
