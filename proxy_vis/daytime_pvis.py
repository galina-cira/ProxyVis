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
import numpy as np
from proxy_vis import solar_zenith_angle

# allowed range of values for visible reflectance
VIS_VALID_MIN  = 0
VIS_VALID_MAX  = 1.3

def vis_disp_sza(
    satellite: str,
    lons: np.ndarray,
    lats: np.ndarray,
    time_info_dt: dt.datetime,
    c02: np.ndarray,
):
    """
    transform Vis channel red (0.64 nm) to make daytime portion of Pvis
    the transformation is done by normalizing Vis by solar zenith angle (SZA) to make
    it visible up to the day/night terminator
    inputs:
        vis: array of vis data, numpy array
        lons: array of lons, numpy array
        lats: array of lats, numpy array
        time_info_dt: datetime object for the time in the middle of the scan
            for example, for Himawary full disk that is start time + 5 minutes
            for 10-min full disk scan

    outputs:
        vis_display : numpy array
                    : adjusted full disk Vis to use with ProxyVis as a daytime part


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