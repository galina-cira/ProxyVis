"""
    The algorithm to normalize full disk GeoProxyVis imagery used in
    the National Weather Service (NWS) Operations

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

# See Chirokova et al. 2023 for details
GAMMA_CORRECTION = 1.0 / 1.5


def normalize_pvis(
    tt_regr_pvis: np.ndarray,
    saved_pvis_min: float,
    saved_pvis_max: float,
    use_saved_params: bool = False,
):
    """
    ProxyVis post-processing to make the image more satble across diferent times
    this normalization is applied to each of the four ProxyVis algorithms


    inputs:
        tt_regr_pvis   :   numpy array
                        :   proxy_vis estimates before normalization

        saved_pvis_min :   float
                        :   saved proxy_vis min values to use

        saved_pvis_max :   float
                        :   saved proxy_vis max values to use

        use_saved_params  : boolean
                          : True: use saved proxy_vis params
                            False: calculate proxy_vis min/max. This will ONLY work for full disk data

    outputs:
        proxy_vis       :   numpy array
                        :   final normalized proxy_vis data

        saved_pvis_min :   float
                        :   proxy_vis min that was used

        saved_pvis_max :   float
                        :   proxy_vis max that was used

    """
    # remove data that are less than zero
    # this is usually just a few points for the GOES-16/GOES_17/Himawari full disk
    less_than_zero_idx = tt_regr_pvis < 0.0
    tt_regr_pvis[less_than_zero_idx] = 0.0

    # normalize ProxyVis values
    if use_saved_params:
        # use pre-saved min/max to normalize
        pvismax = saved_pvis_max
        pvismin = saved_pvis_min
    else:
        pvismax = np.nanmax(tt_regr_pvis)
        pvismin = np.nanmin(tt_regr_pvis)

    regr_pvis_norm = (tt_regr_pvis - pvismin) / (pvismax - pvismin)

    # Remove negative and NaN values. This is required to remove the black dots
    # that could overvise be seen, especially in Himawari and Meteosat second generations imagery
    # The negative values that are being removed are related to the noise in 3.9 um channel
    regr_pvis_norm[less_than_zero_idx] = 0.0
    regr_pvis_norm[np.isnan(regr_pvis_norm)] = pvismax

    # estimate the final ProxyVis product. This step is needed to make
    # nighttime and daytime parts of the GeoProxyVis imagery look more similar.
    proxy_vis = regr_pvis_norm**GAMMA_CORRECTION

    return proxy_vis, pvismin, pvismax
