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

import numpy as np

from proxy_vis import norm_pvis, saved_pvis_min_max

# all clouds
# non-normalized regression parameters
PARAMS = np.array([-3.04491517e-13, 8.16054268e-01])


def calculate_pvis_simple_one_eq(
    satellite: str, c07: np.ndarray, use_saved_params: bool = False
):
    """
    Single-channel one-regression ProxyVis algorithm


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
