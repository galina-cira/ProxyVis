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
):
    """
    This is NOT the main algorithm, see main2 for the main
    multi-channel multiple-regression operational ProxyVis algorithm


    inputs:
        satellite               : string
                                : satellite to process

        all_inp_bands_data_dict :   dictionary
                                :   dictionary is used to allow common interface
                                    for different versions that use different channels
                                    keys - GOES/HW channels ('band02', etc),
                                    values - corresponding data arrays

        use_saved_params        : set True to use saved values; set False to estimate min/max from data


    outputs:

        proxy_vis               : numpy array
                                : full disk nighttime ProxyVis at original 2km resolution; numpy array

        tt_for_regr             : numpy array
                                : intermediate values only used for development

        pvismin                : float
                                : ProxyVis min estimated from data

        pvismax                : float
                                : ProxyVis min estimated from data


    """

    saved_pvis_min, saved_pvis_max = saved_pvis_min_max.lookup_range(satellite)

    # calculate ProxyVis from individual bands

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
