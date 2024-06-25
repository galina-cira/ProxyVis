"""
    Algorithm: multi-channel single-regressions ProxyVis algorithm; this is NOT
    the main operational algorithm see calculate_pvis_main_two_eq for the main
    operational algorithm

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

from typing import Tuple

import numpy as np

from proxy_vis import norm_pvis, saved_pvis_min_max

# Min value to use with log function to avoid large negative values
MIN_LOG_VAL = 3e-05

# non-normalized regressionparameters
PARAMS = np.array([-3.44571376e-02, -2.50124844e-02, -2.95821592e-13, 8.82291378e-01])


def calculate_pvis_main_one_eq(
    satellite: str,
    c07: np.ndarray,
    c11: np.ndarray,
    c13: np.ndarray,
    c15: np.ndarray,
    use_saved_params: bool = False,
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """ProxyVis Multi-Channel Single-Regression Equation

    This is NOT the main operational algorithm, see nighttime_pvis_main_two_eq
        for the main multi-channel multiple-regression operational ProxyVis algorithm

    Args:
        satellite (str): Name of the satellite ("goes16", "goes17", "goes18",
            "himawari8", "himawari9", "meteosat-9", "meteosat-11")
        c07 (np.ndarray): IR channel similar to GOES-16 C07 channel data
        c11 (np.ndarray): IR channel similar to GOES-16 C11 channel data
        c13 (np.ndarray): IR channel similar to GOES-16 C13 channel data
        c15 (np.ndarray): IR channel similar to GOES-16 C15 channel data
        use_saved_params (bool):  Whether to use static/dynamic normalization
            for ProxyVis. Dynamic normalization works best when estimated from the full
            disk GEO data. For sub-sector ProxyVis images it is recommended to use saved
            normalization parameters.

    Returns:
        proxy_vis (np.ndarray): Full disk nighttime ProxyVis at the original 
            2-km resolution.
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
    tt2 = np.full_like(c07, np.nan)
    tt3 = np.full_like(c07, np.nan)
    tt = np.full_like(c07, np.nan)

    # define the min value to use in log function
    diff1 = abs(c11 - c07)
    diff1[diff1 < MIN_LOG_VAL] = MIN_LOG_VAL
    tt0 = PARAMS[-1]
    tt1 = PARAMS[2] * (c07**5)
    tt2 = PARAMS[1] * (np.log(diff1))
    tt3 = PARAMS[0] * ((abs(c13 - c15)) ** 0.4)

    # assemble full array
    tt = tt0 + tt1 + tt2 + tt3
    regr_pvis = np.copy(tt)
    proxy_vis, pvismin, pvismax = norm_pvis.normalize_pvis(
        tt, saved_pvis_min, saved_pvis_max, use_saved_params
    )

    return proxy_vis, regr_pvis, pvismin, pvismax
