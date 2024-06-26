"""
    Vis algorithm:the main algorithm to create adjusted Vis daytime imagery
    used to generate GeoProxyVis in the National Weather Service (NWS) Operations

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
"""

import datetime as dt
from typing import Tuple

import numpy as np

from proxy_vis import solar_zenith_angle

# allowed range of values for visible reflectance
VIS_VALID_MIN = 0
VIS_VALID_MAX = 1.3


def vis_disp_sza(
    satellite: str,
    lons: np.ndarray,
    lats: np.ndarray,
    time_info_dt: dt.datetime,
    c02: np.ndarray,
) -> Tuple[np.ndarray, float, float]:
    """Generate adjusted Vis data.

    Transform Vis channel (Red, 0.64 nm, Ch 02 in GOES-16) to make daytime
    portion of GeoProxyVis.  The transformation is done by normalizing Vis by solar
    zenith angle (SZA) to make it visible up to the day/night terminator.

    Args:
        satellite (str): Name of the satellite ("goes16", "goes17", "goes18",
            "himawari8", "himawari9", "meteosat-9", "meteosat-11")
        lons (np.ndarray): array of lons
        lats (np.ndarray): array of lats
        time_info_dt (dt.datetime): Time of the MIDDLE of the scan. This is NOT
            the GEO satellite data timestamp. The GEO timestamp is at the start
            of the scan. This is the time of the middle of the scan to calculate SZA most
            relevant for the full disk. For example, for GOES-16/18, and Himawari8/9 
            10-min full disk scan this is teh timestamp in the filename +5 minutes. 
        c02 (np.ndarray): Visible channel similar to GOES-16 C02 channel data

    Returns:
        vis_disp (np.ndarray): Adjusted full disk Vis to use for the daytime
            portion of GeoProxyVis 
        vismin (float): The minimum value used to normalize the
            Vis data.  
        vismax (float): The maximum value used to normalize the Vis data.
    """

    # get arrey of solar zenith angle (SZA) values for the full disk and Vis resolution
    # usually that will be 0.5 km grid
    rsun_zen = solar_zenith_angle.sza_pysolar(time_info_dt, lons, lats)
    # convert to radians
    rsun_zen_pos_rad = np.deg2rad(rsun_zen)

    # ensure visible data are between 0 and 1.3
    too_big = (np.isfinite(c02)) & (c02 > 1.3)
    too_small = (np.isfinite(c02)) & (c02 < 0)
    c02[too_big] = VIS_VALID_MAX
    c02[too_small] = VIS_VALID_MIN

    # calculate Vis values that will be combined with the nighttime part
    vis_disp = c02 / np.cos(rsun_zen_pos_rad)
    vis_disp **= 0.5

    vismin = np.nanmin(c02)
    vismax = np.nanmax(c02)

    return vis_disp, vismin, vismax
