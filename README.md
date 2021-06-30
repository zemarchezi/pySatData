# Basic documentation for pySatData

## Requirements


Python 3.8+

Pandas, NumPy, wget, requests, cdflib, pytplot

Installation: ```pip instal -e .```

***
## config_file.json

```pysatdata/resources/config_file.json```

This file sets the http directory for downloading data and the local directory to download data.

The local directory files are organized as http directory, i.e.:
* RBSP: HOME/sat_data/rbsp/rbspa/l2/ect/rept/sectors/rel03/YYYY/filename.cdf

This files handles with the different http subpaths for the different levels and intruments in each probe.

If some change is needed, it can be done in this file.
***

```python
import json
from pysatdata.loaders.load import *

trange=['2021-05-26', '2021-05-30']

varss_rept = load_sat(trange=trange, satellite='rbsp',
                     probe=['a'], level='2', rel='rel03',
                     instrument='rept',datatype='sectors',
                     config_file='./pysatdata/resources/config_file.json', downloadonly=False,
                     usePandas=False, usePyTplot=True)
```