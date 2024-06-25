"""
    ProxyVis algoritm:simple single regression

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
    #############################################################################
"""

from typing import Tuple

import numpy as np

from proxy_vis import norm_pvis, saved_pvis_min_max

# all clouds
# non-normalized regression parameters
PARAMS = np.array([-3.04491517e-13, 8.16054268e-01])


def calculate_pvis_simple_one_eq(
    satellite: str, c07: np.ndarray, use_saved_params: bool = False
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """ProxyVis Single Channel One Regression Algorithm

    Args:
        satellite (str): Name of the satellite ("goes16", "goes17", "goes18",
            "himawari8", "himawari9", "meteosat-9", "meteosat-11")
        c07 (np.ndarray): IR channel similar to GOES-16 C07 channel data
        use_saved_params (bool):  Whether to use static/dynamic normalization for
            for ProxyVis. Dynamic normalization works best when estimated from the full
            disk GEO data. For sub-sector ProxyVis images it is recommended to use saved
            normalization parameters.

    Returns:
        proxy_vis (np.ndarray): Full disk nighttime ProxyVis at the original
            2km resolution.
        regr_pvis (np.ndarray): Intermediate values only used for development.
        pvismin (float): The minimum value used to normalize the ProxyVis data.
            If use_saved_params was True, then this value was looked up based
            on the satellite.  Otherwise, this is calculated based on the data
            in data_dict.
        pvismax (float): The maximum value used to normalize the ProxyVis data.
            If use_saved_params was True, then this value was looked up based
            on the satellite.  Otherwise, this is calculated based on the data
            in data_dict.
    """

    saved_pvis_min, saved_pvis_max = saved_pvis_min_max.lookup_range(satellite)

    # calculate ProxyVis from individual channels

    # define arrays
    tt0 = np.full_like(c07, np.nan)
    tt1 = np.full_like(c07, np.nan)
    tt = np.full_like(c07, np.nan)

    # all clouds
    tt0 = PARAMS[-1]
    tt1 = PARAMS[0] * (c07**5)

    # assemble full array
    tt = tt0 + tt1
    regr_pvis = np.copy(tt)
    proxy_vis, pvismin, pvismax = norm_pvis.normalize_pvis(
        tt, saved_pvis_min, saved_pvis_max, use_saved_params
    )

    return proxy_vis, regr_pvis, pvismin, pvismax
