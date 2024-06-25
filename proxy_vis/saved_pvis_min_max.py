"""
    This module provides the saved min and max values for ProxyVis for all
    satellites described in Chirokova et al. (2023). Originally ProxyVis was
    normalized using full-disk min and max values. This created various issues
    with implementing ProxyVis in operations. Thus, we estimated typical min and
    max values that can be used instead of dynamic normalization values.

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

SAVED_PVIS_MIN = {}
SAVED_PVIS_MAX = {}

SAVED_PVIS_MIN["goes16"] = 0.0
SAVED_PVIS_MAX["goes16"] = 0.78

SAVED_PVIS_MIN["goes17"] = 0.0
SAVED_PVIS_MAX["goes17"] = 0.84

SAVED_PVIS_MIN["goes18"] = 0.0
SAVED_PVIS_MAX["goes18"] = 0.84

SAVED_PVIS_MIN["himawari8"] = 0.0
SAVED_PVIS_MAX["himawari8"] = 0.79

SAVED_PVIS_MIN["himawari9"] = 0.0
SAVED_PVIS_MAX["himawari9"] = 0.79

SAVED_PVIS_MIN["meteosat-9"] = 0.0
SAVED_PVIS_MAX["meteosat-9"] = 0.78

SAVED_PVIS_MIN["meteosat-10"] = 0.0
SAVED_PVIS_MAX["meteosat-10"] = 0.78

SAVED_PVIS_MIN["meteosat-11"] = 0.0
SAVED_PVIS_MAX["meteosat-11"] = 0.78


def lookup_range(satellite: str) -> Tuple[float, float]:
    """Looks up the min/max normalization value for the given satellite.

    Args:
        satellite (str): The name of the satellite. Must be one of: goes16, goes17
            goes18, himawari8, himawari9, meteosat-9, meteosat-10, or meteosat-11.

    Raises:
        KeyError: If the given satellite is not in the list of supported satellites,
            then a KeyError will be raised.

    Returns:
        saved_pvis_min(float): The looked up min.
        saved_pvis_max(float): The looked up max.
    """
    sanitized = satellite.lower().strip()

    try:
        saved_pvis_min = SAVED_PVIS_MIN[sanitized]
        saved_pvis_max = SAVED_PVIS_MAX[sanitized]
    except KeyError as exc:
        raise KeyError(
            f"Creating ProxyVis for satellite {satellite} is not implemented."
        ) from exc

    return saved_pvis_min, saved_pvis_max
