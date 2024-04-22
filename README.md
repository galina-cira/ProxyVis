# ProxyVis

## High Level API Usage

This is a high-level API for generating GeoProxyVis imagery. The GeoProxyVis imagery
combines a nighttime infrared (IR)-based proxy for Visible (VIS) imagery at nighttime and adjusted VIS imagery at
daytime to create combined day-night GeoProxyVis imagery. The algorithm is described in

Chirokova, G., J. A. Knaff, M. J. Brennan, R. T. DeMaria, M. Bozeman, S. N.
Stevenson, J. L. Beven, E. S. Blake, A. Brammer, J. W. Darlow, M. DeMaria, S.
D. Miller, C. J. Slocum, Debra Molenar, and D. W. Hillger, 2023: ProxyVis—A
Proxy for Nighttime Visible Imagery Applicable to Geostationary Satellite
Observations. Wea. Forecasting, 38, 2527–2550,
https://doi.org/10.1175/WAF-D-23-0038.1.

A couple of important notes:
1)  The ProxyVis algorithm is only valid during the nighttime (defined as solar
zenith angle (SZA) > 89 deg. Most of the signal in ProxyVis is from 3.9 um
channel. Since 3.9 um channel has both reflective and emissive components, it
behaves differently during daytime and nighttime. The ProxyVis algorithm was
developed for nighttime and has not been verified or tested in any way for
daytime values. Thus, the ProxyVis values are undefined during daytime, for
solar zenith angles < 89.0 deg 
2) ProxyVis is generated at 2 km, the native resolution for IR channles.
However, in order to create a combined full disk image without degrading the
resolution of Vis channel, for operational use the 2 km ProxyVis is regrided to
0.5 km Vis channel resolution.  This software provides full disk combined day/night
GeoProxyVis imagery at both 0.5 km and 2 km resolution. 

Description:
This high level API will take a dictionary of IR and visible input
geostationary satellite data, a string representing the ProxyVis function to
use, and a mapping from the channels in the input dictionary to the function
arguments and will produce a composite ProxyVis image.

The input dictionary should take the form `data_dict[field_type][channel_name]`
where `channel_name` is the name of the channel (such as `C07` for channel 7
GOES data) and `field_type` is either `bt_temp` if the data is IR data or
`radiances` if the data is visible data.

Constants that map each supported satellite's channel names in the `data_dict`
to each ProxyVis function's arguments is provided by the `channel_to_arg_maps`
module for convenience.  The provided dictionaries assume the channel names
match those provided by Satpy. In your own code, the dictionary mapping channel
names to function arguments should probably come from a config file.

Example:
```python
from proxy_vis import channel_to_arg_maps, combine_dn_pvis
# Read GOES data into a dictionary that looks like:
# data_dict["C02"]["radiances"]['C02']
# data_dict["bt_temp"]["C07"]
# data_dict["bt_temp"]["C11"]
# data_dict["bt_temp"]["C13"]
# data_dict["bt_temp"]["C15"]
data_dict, data_time, vis_lons, vis_lats, ir_lons, ir_lats = your_reader(goes_filenames)

# Most GOES fulldisk imagery is produced every 10 minutes.
# This is used by the lower level code to generate a day/night mask at the 
# center time of the scan.
GOES_MINUTE_INTERVAL = 10
use_saved_params = True

# Get the composite ProxyVis data as well as the min/max value used during
# normalization of the ProxyVis data.
pvis_combined_05km, pvis_combined_2km, pvismin, pvismax = combine_dn_pvis.get_all_vis_pvis(
    "goes16", # Name of the satellite ("goes16", "goes17", "goes18", "himawari8", "himawari9", "meteosat-9", "meteosat-11")
    data_time, # Time of the provided data, provided as python datetime object
    data_dict, # Dictionary of IR and visible data
    pvis_channel_to_arg_maps.ABI_MAIN, # Dict mapping data_dict entries to pvis function args
    dvis_channel_to_arg_maps.ABI_MAIN, # Dict mapping data_dict entries to daytime vis function args
    ir_lons, # Longitudes at the IR resolution (2km)
    ir_lats, # Latitudes at the IR resolution (2km)
    vis_lons, # Longitudes at the visible resolution (0.5km)
    vis_lats, # Latitudes at the visible resolution (0.5km)
    SATELLITE_MINUTE_INTERVAL, # Length in minutes of the scan of the full disk
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

##########################################################################

This code is part of the ProxyVis processing written by:
Galina.Chirokova@colostate.edu; Robert.DeMaria@colostate.edu,
Alan Brammer

Copyright (C) 2024  Galina Chirokova, Robert DeMaria, Alan Brammer

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

