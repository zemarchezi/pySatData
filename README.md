# Basic documentation for pySatData


## Requirements


Python 3.8+

1. Clone or download the pySatData
2. From the pySatData directory: 

Installation: ```pip instal -e .```

***
## config_file.json

If some change in downloding directories is needed, it can be done in this file.

```pysatdata/resources/config_file.json```

This file sets the http directory for downloading data and the local directory to save the downloaded data.

The local directory files are organized as http directory, i.e.:
* RBSP: HOME/sat_data/rbsp/rbspa/l2/ect/rept/sectors/rel03/YYYY/filename.cdf

This files handles with the different http subpaths for the different levels and intruments in each probe.

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


See ```plot_interpFlux_RBSP.py``` for an example on plotting the interpolated electron flux 
