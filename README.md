# ProxyVis and GeoProxyVis Imagery for Geostationary Satellite Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12557264.svg)](https://doi.org/10.5281/zenodo.12560629)

This package provides high-level and low-level APIs for generating GeoProxyVis
imagery. https://rammb2.cira.colostate.edu/research/goes-r-research/proxyvis/ 

GeoProxyVis imagery combines ProxyVis - a nighttime infrared (IR)-based proxy
for Visible (VIS, 0.64 um) imagery at nighttime, and adjusted VIS imagery at
daytime to create a combined day-night GeoProxyVis imagery. The algorithm is
described in:

Chirokova, G., J. A. Knaff, M. J. Brennan, R. T. DeMaria, M. Bozeman, S. N.
Stevenson, J. L. Beven, E. S. Blake, A. Brammer, J. W. Darlow, M. DeMaria, S.
D. Miller, C. J. Slocum, D. Molenar, and D. W. Hillger, 2023: ProxyVis—A
Proxy for Nighttime Visible Imagery Applicable to Geostationary Satellite
Observations. Wea. Forecasting, 38, 2527–2550,
https://doi.org/10.1175/WAF-D-23-0038.1.


A couple of important notes:
1)  The ProxyVis algorithm is only valid during the nighttime (defined as solar
zenith angle (SZA) > 89 deg). Most of the signal in ProxyVis is from the 3.9-um
channel. Since the 3.9-um channel has both reflective and emissive components,
it behaves differently during daytime and nighttime. The ProxyVis algorithm was
developed for nighttime and has not been verified or tested in any way for
daytime values. Thus, the ProxyVis values are undefined during daytime, for
solar zenith angles < 89.0 deg 
2) ProxyVis is generated at 2 km, the native resolution for IR channles.
However, in order to create a combined full disk image we need both the VIS and
ProxyVis parts of the full disk image at the same resolution. To achieve that,
we need to either regrid VIS to 2-km or regrid ProxyVis to 0.5-km resolution.
For the National Weather Service (NWS) operational use, ProxyVis is regridded
to 0.5-km so that GeoProxyVis can provide VIS data at the native 0.5-km
resolution. The high-level API can produce combined full disk day/night
GeoProxyVis imagery at both 0.5 km and 2 km resolutions. 

## High Level API Usage 
The high level API will take a dictionary of IR and VIS input geostationary
satellite data, a string representing the ProxyVis function to use, a string
representing the VIS function to use, a dictionary mapping channel names to
ProxyVis function arguments, and a dictionary mapping channel names to VIS
function arguments. This will produce a full disk day/night GeoProxyVis image.

The input dictionary should take the form `data_dict[field_type][channel_name]`
where `channel_name` is the name of the channel (such as `C07` for channel 7
GOES data) and `field_type` is either `bt_temp` if the data is IR data or
`radiances` if the data is VIS data.

The available ProxyVis function names are one of 
    `nighttime_pvis_main_two_eq`,
    `nighttime_pvis_main_one_eq`, 
    `nighttime_pvis_simple_two_eq`,
    `nighttime_pvis_simple_one_eq`.  
See the modules with the same names to examine the ProxyVis functions. The
`nighttime_pvis_main_two_eq` is the main ProxyVis algorithm used in NWS
operations. 

There is currently only one VIS function provided: `vis_disp_sza`.

Constants that map each supported satellite's channel names in the `data_dict`
to each ProxyVis/VIS  function's arguments are provided by the
`channel_to_arg_maps` module for convenience.  The provided dictionaries assume
the channel names match those provided by Satpy. In your own code, the
dictionary mapping channel names to function arguments should probably come
from a config file.

Example:
```python
from proxy_vis import channel_to_arg_maps, combine_dn_pvis
# Read GOES data into a dictionary that looks like:
# data_dict["radiances"]["C02"]
# data_dict["bt_temp"]["C07"]
# data_dict["bt_temp"]["C11"]
# data_dict["bt_temp"]["C13"]
# data_dict["bt_temp"]["C15"]
data_dict, data_time, vis_lons, vis_lats, ir_lons, ir_lats = your_reader(goes_filenames)

# Current GOES-16 and GOES-18 full disk imagery is produced every 10 minutes.
# This is used by the lower level code to generate a day/night mask at the 
# center time of the scan.
GOES_MINUTE_INTERVAL = 10

# Use a saved set of min/max values for each satellite. If this is set
# to false, the min/max values will be computed from the full disk image data at runtime.
use_saved_params = True

# Get the full disk GeoProxyVis data as well as the min/max values used to
# normalize the ProxyVis data.
pvis_combined_05km, pvis_combined_2km, pvismin, pvismax = combine_dn_pvis.get_all_vis_pvis(
    "goes16", # Name of the satellite ("goes16", "goes17", "goes18", "himawari8", "himawari9", "meteosat-9", "meteosat-11")
    data_time, # Time of the provided data, given as python datetime object
    data_dict, # Dictionary of IR and VIS data
    channel_to_arg_maps.ABI_MAIN, # Dict mapping data_dict entries to pvis function args
    channel_to_arg_maps.ABI_VIS, # Dict mapping data_dict entries to daytime vis function args
    ir_lons, # Longitudes at IR resolution (2km)
    ir_lats, # Latitudes at IR resolution (2km)
    vis_lons, # Longitudes at VIS resolution (0.5km)
    vis_lats, # Latitudes at VIS resolution (0.5km)
    GOES_MINUTE_INTERVAL, # Length in minutes of the scan of the full disk
    "nighttime_pvis_main_two_eq", # Desired ProxyVis function name
    "vis_disp_sza", # Desired daytime Vis function name
    use_saved_params=use_saved_params, # Whether to use saved/dynamic normalization for ProxyVis
    output_resolution="both", # Must be one of "2.0km", "0.5km", or "both". If "2.0km" or "0.5km"
        # are selected then either the returned 2km or 0.5km array will be set to None.
) 
```

## Himawari Example Script
An example script that makes use of the high level API can be found in the
`examples/ahi_example.py` file. This example assumes you have the `proxy_vis`,
`satpy`, `pyresample`, and `matplotlib` libraries installed in your environment.
To run, place `B03`, `B07`, `B11`, `B13`, and `B15` Himawari-9 data for a single
time in the `examples/input/` directory. 

To run the code, you can use the command:
`python examples/ahi_example.py`

This will produce a plot of ProxyVis data in the `examples/output/` directory.

## License
##########################################################################

This code is part of the ProxyVis processing written by:
Galina.Chirokova@colostate.edu; Robert.DeMaria@colostate.edu,
Alan Brammer

Copyright (C) 2018 - 2024  Galina Chirokova, Robert DeMaria, Alan Brammer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
##########################################################################

