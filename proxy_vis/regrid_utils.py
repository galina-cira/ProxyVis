"""
    Utilities to regrid data between 2 km and 0.5 km grids

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

from typing import Optional

import numpy as np
import pyresample as pr

# pyresample parameters
EPSILON = 0.5
DEFAULT_ROI_M = 10000


def regrid_nn(
    src_lons: np.ndarray,
    src_lats: np.ndarray,
    src_data: np.ndarray,
    dst_lons: np.ndarray,
    dst_lats: np.ndarray,
    radius_of_influence_m: int = DEFAULT_ROI_M,
    fill_value: Optional[float] = None,
    epsilon: float = 0.5,
    should_fill_mask: bool = True,
) -> np.ndarray:
    """Regrids data using nearest-neighbor interpolation.

    This is a convenience wrapper for pyresample's kdtree nearest neighbor
    interpolation.  By default, a missing value of NaN will be used.

    Args:
        src_lons (np.ndarray): Array of longitudes for the src_data.
        src_lats (np.ndarray): Array of latitudes for the src_data.
        src_data (np.ndarray): Source data to regrid.
        dst_lons (np.ndarray): Array of longitudes to regrid to.
        dst_lats (np.ndarray): Array of latitudes to regrid to.
        radius_of_influence_m (int, optional): Radius of influence in meters to
            use for nearest neighbor interpolation.  Only source points within
            this radius from the destination points are considered. Defaults to
            10000.
        fill_value (float, optional): The fill value to use for
            missing points.  If None, NaN will be used for missing points.
            Defaults to None.
        epsilon (float, optional): Threshold of error allowed in nearest
            neighbor distance calculation, this can greatly speed up the
            interpolation calculation. See pyresample documentation for details.
            Defaults to 0.5.
        should_fill_mask (bool, optional): Pyresample returns masked arrays, but
            this can cause problems for some routines that don't support masked
            arrays (such as imshow).  If this parameter is true, then the masked
            array will be converted to a normal array where missing data is set
            to the fill_value(NaN is used if the fill_value is None). Defaults
            to True.

    Returns:
        result (np.ndarray): The data regridded to the destination grid.
    """
    src_area = pr.geometry.SwathDefinition(lons=src_lons, lats=src_lats)
    dst_area = pr.geometry.SwathDefinition(lons=dst_lons, lats=dst_lats)

    result = pr.kd_tree.resample_nearest(
        src_area,
        src_data,
        dst_area,
        radius_of_influence=radius_of_influence_m,
        fill_value=fill_value,
        epsilon=epsilon,
    )

    if should_fill_mask:
        fill_to_use = fill_value
        if fill_to_use is None:
            fill_to_use = np.nan
        result = result.filled(fill_to_use)

    return result
