## Welcome to ```STATION_ASF```, a NEMO setup which mimics an off-shore station dedicated to the measure of Air-Sea Fluxes (ASF).



### Important namelist parameters

* ```rn_dept1@namusr_def:``` depth (m) at which the prescribed SST is taken (i.e. depth of first T-point); important due to impact on warm-layer estimate, the deeper, the more pronounced!

* ```rn_lat1d,rn_lon1d@namc1d:``` fixed coordinates of the location of the station (buoy, platform, etc).

* ```namsbc_blk:``` to be filled carefully, just as for "C1D", the prescribed surface ATMOSPHERIC state (files) are time series of shape 3x3 in space

* ```namsbc_sas:``` to be filled carefully, just as for "C1D", the prescribed surface OCEAN state (files) are time series of shape 3x3 in space
