"""
    Utilities to create day/night mask for given time and lat/lon arrays

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
from pysolar import solar


def sza_pysolar(
    sza_datetime: dt.datetime,
    lons_e: np.ndarray,
    lats_n: np.ndarray,
    tzinfo: dt.timezone = dt.timezone.utc,
):
    """
    calculate solar zenith angle array using pysolar
    inputs:
        sza_datatime: datetime object
        lon_e: longitutes East (numpy array)
        lat_n: latitides north (numpy array)
        tzinfo  : datetime time zone info (optional)
    outputs:
        sza: array of sza (numpy array)
    """

    # calculate full disk solar zenith angle array
    sza = 90.0 - solar.get_altitude(lats_n, lons_e, sza_datetime.replace(tzinfo=tzinfo))

    return sza


def sza_to_mask(sza, mask_night=True, sza_thresh=90.0):
    """
    make overlapping sza mask with sza = threshold assigned to both day and night
    inputs:
        sza: array of solar zenith angles (numpy array)
        mask_night: True/False - set night to NaN is True
                               - set day to NaN if False
        sza_thresh: default threshold to use for SZA
            Note: many GEO day/night applications require sza_thresh = 89.0
            to avoid dividing by zero
    outputs:
        sza: array of sza with day or night set to NaN
    """

    if mask_night:
        sza[sza > sza_thresh] = np.nan
        sza[sza <= sza_thresh] = 1.0
    else:
        sza[sza <= sza_thresh] = np.nan
        sza[sza > sza_thresh] = 1.0

    return sza
