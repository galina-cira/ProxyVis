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
) -> np.ndarray:
    """Calculates solar zenith angle array using pysolar

    Args:
        sza_datetime (dt.datetime): Time to generate SZA for
        lon_e (np.ndarray): longitutes East
        lat_n (np.ndarray): latitides North
        tzinfo (dt.timezone): Time zone to use for sza_datetime. Default is UTC

    Returns:
        sza (np.ndarray): Solar zenith angle for the given time
    """

    # calculate full disk solar zenith angle array
    sza = 90.0 - solar.get_altitude(lats_n, lons_e, sza_datetime.replace(tzinfo=tzinfo))

    return sza


def sza_to_mask(sza: np.ndarray, mask_night: bool = True, sza_thresh: float = 90.0):
    """Generates a mask from the given solar zenith angles

    Makes overlapping sza mask with sza = threshold assigned to both day and
    night data.

    Args:
        sza (np.ndarray): Array of solar zenith angles
        mask_night (bool): Set night to NaN if True, otherwise set day to NaN.
            Default is True
        sza_thresh (float): Threshold to use for SZA. Many GEO day/night
            applications require sza_thresh = 89.0 deg to avoid dividing by zero.
            Default is 90.0

    Returns:
        sza (np.ndarray): Array of SZA with either day or night set to NaN
    """

    if mask_night:
        sza[sza > sza_thresh] = np.nan
        sza[sza <= sza_thresh] = 1.0
    else:
        sza[sza <= sza_thresh] = np.nan
        sza[sza > sza_thresh] = 1.0

    return sza
