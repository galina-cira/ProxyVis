"""
    This is the main module to esimate both daytime Vis and nighttime ProxyVis
    parts of the GeoProxyVis imagery. Your code shoud call the
    get_all_vis_pvis function.

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
from typing import Dict, Tuple

import numpy as np

from proxy_vis import (
    daytime_pvis,
    nighttime_pvis_main_one_eq,
    nighttime_pvis_main_two_eq,
    nighttime_pvis_simple_one_eq,
    nighttime_pvis_simple_two_eq,
    solar_zenith_angle,
)
from proxy_vis.regrid_utils import regrid_nn

# Input data format
DataDict = Dict[str, Dict[str, np.ndarray]]
"""Dictionary whose keys are channel names and values are dictionaries of data.

Dictionary of data read from a file.  The keys should be channel names and the
values should be a sub-dictionary of calibrated data.  The keys for the 
sub-dictionary should be the type of calibration: one of "bt_temp", "radiances".
"""

# Available ProxyVis functions
PVIS_FUNC_LOOKUP = {
    "nighttime_pvis_main_two_eq": nighttime_pvis_main_two_eq.calculate_pvis_main_two_eq,
    "nighttime_pvis_main_one_eq": nighttime_pvis_main_one_eq.calculate_pvis_main_one_eq,
    "nighttime_pvis_simple_two_eq": nighttime_pvis_simple_two_eq.calculate_pvis_simple_two_eq,
    "nighttime_pvis_simple_one_eq": nighttime_pvis_simple_one_eq.calculate_pvis_simple_one_eq,
}

# Available Vis functions
VIS_FUNC_LOOKUP = {"vis_disp_sza": daytime_pvis.vis_disp_sza}

# supported output resolutions
OUTPUT_RES_2KM = "2.0km"
OUTPUT_RES_05KM = "0.5km"
OUTPUT_RES_BOTH = "both"
VALID_OUTPUT_RES_KEYWORDS = (OUTPUT_RES_2KM, OUTPUT_RES_05KM, OUTPUT_RES_BOTH)

# threshold to use to separate daytime and nighttime parts
SZA_THRESHOLD = 89.0
# radiances in Vis channel are between 0 and 1.3
# the ProxyVis has the values from 0 to 1. We scale that by 1.3 to
# get the same range of values before merging Vis and ProxyVis parts
VIS_SCALING_FACTOR = 1.3


def get_all_vis_pvis(
    satellite: str,
    time_dt: dt.datetime,
    data_dict: DataDict,
    pvis_data_to_args_map: Dict[str, str],
    dvis_data_to_args_map: Dict[str, str],
    geo_lons_2km: np.ndarray,
    geo_lats_2km: np.ndarray,
    geo_lons_05km: np.ndarray,
    geo_lats_05km: np.ndarray,
    minutes_interval: int,
    pvis_alg: str,
    dvis_alg: str,
    use_saved_params: bool,
    output_resolution: str = OUTPUT_RES_BOTH,
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """Generate combined day/night full disk GeoProxyVis imagery.

    Args:
        satellite (str): Name of the satellite ("goes16", "goes17", "goes18",
            "himawari8", "himawari9", "meteosat-9", "meteosat-11")
        time_dt (dt.datetime): Timestamp of the input satellite data.
        data_dict (DataDict): Dictionary of IR and Vis data.
        pvis_data_to_args_map (Dict[str, str]): Dict mapping data_dict entries
            to pvis function args
        dvis_data_to_args_map (Dict[str, str]): Dict mapping data_dict entries
            to daytime vis function args
        geo_lons_2km (np.ndarray): GEO satellite longitudes at IR resolution (2km)
        geo_lats_2km (np.ndarray): GEO satellite latitudes at IR resolution (2km)
        geo_lons_05km (np.ndarray): GEO satellite longitudes at Vis resolution (0.5km)
        geo_lats_05km (np.ndarray): GEO satellite latitudes at Vis resolution (0.5km)
        minutes_interval (int): Length in minutes of the scan of the full disk
        pvis_alg (str): User requested ProxyVis function name:
            ("nighttime_pvis_main_two_eq", "nighttime_pvis_main_one_eq",
            "nighttime_pvis_simple_two_eq", "nighttime_pvis_simple_one_eq")
        dvis_alg (str): User requested daytime Vis function name, currently only
            "vis_disp_sza" is supported.
        use_saved_params (bool): Whether to use static/dynamic normalization
            for ProxyVis. Dynamic normalization works best when estimated from the full
            disk GEO data. For sub-sector ProxyVis images it is recommended to use saved
            normalization parameters.
        output_resolution (str, optional): Must be one of "2.0km", "0.5km", or
            "both". If "2.0km" or "0.5km" are selected then either the returned
            2km or 0.5km array will be set to None. Defaults to "both".

    Returns:
        pvis_combined_05km (np.ndarray): Array containing the combined day/night
            GeoProxyVis image at 0.5km resolution.  If output_resolution was
            set to "2.0km" then this will be None.
        pvis_combined_2km (np.ndarray): Array containing the combined day/night
            GeoProxyVis image at 2km resolution.  If output_resolution was
            set to "0.5km" then this will be None.
        pvismin (float): The minimum value used to normalize the ProxyVis data.
            If use_saved_params was True, then this value was looked up based
            on the satellite.  Otherwise, this is calculated based on the data
            in data_dict.
        pvismax (float): The maximum value used to normalize the ProxyVis data.
            If use_saved_params was True, then this value was looked up based
            on the satellite.  Otherwise, this is calculated based on the data
            in data_dict.
    """
    sanitized_out_res = output_resolution.lower().strip()
    _validate_output_res_keyword(sanitized_out_res)

    # get datetime for the middle of the scan. This is used to
    # as a timestamp to claulate SZA and combine day- and night- time parts of the image
    # minutes_interval == 10 min for current Himawari, GOES16/17 (2022), and
    # 15 minutes for GOES16/17 prior to 2019
    dt_time_middle_of_scan = time_dt + dt.timedelta(minutes=minutes_interval / 2.0)

    # get 2 km nighttime proxyvis for the full disk
    # NOTE: invalid values will be produced for daytime and removed while combining
    # day- and night-time parts
    pvis_2km, pvismin, pvismax = calculate_pvis(
        satellite,
        data_dict,
        pvis_data_to_args_map,
        pvis_alg,
        use_saved_params,
    )

    # create adjusted Vis data array
    vis_disp_05km, _vismin, _vismax = calculate_vis(
        satellite,
        data_dict,
        dvis_data_to_args_map,
        dvis_alg,
        geo_lons_05km,
        geo_lats_05km,
        dt_time_middle_of_scan,
    )

    pvis_combined_05km = None
    pvis_combined_2km = None

    if sanitized_out_res in (OUTPUT_RES_05KM, OUTPUT_RES_BOTH):
        # get 05 km pvis
        pvis_05km = regrid_nn(
            geo_lons_2km, geo_lats_2km, pvis_2km, geo_lons_05km, geo_lats_05km
        )

        # get combined dn image at 0.5 km
        pvis_combined_05km = apply_dn_mask(
            pvis_05km,
            vis_disp_05km,
            dt_time_middle_of_scan,
            geo_lats_05km,
            geo_lons_05km,
        )

    if sanitized_out_res in (OUTPUT_RES_2KM, OUTPUT_RES_BOTH):
        # get 2.0 km dvis
        vis_disp_2km = regrid_nn(
            geo_lons_05km, geo_lats_05km, vis_disp_05km, geo_lons_2km, geo_lats_2km
        )

        # get combined dn image at 2.0 km
        pvis_combined_2km = apply_dn_mask(
            pvis_2km, vis_disp_2km, dt_time_middle_of_scan, geo_lats_2km, geo_lons_2km
        )

    return pvis_combined_05km, pvis_combined_2km, pvismin, pvismax


def _validate_output_res_keyword(sanitized_output_res: str) -> None:
    if sanitized_output_res not in VALID_OUTPUT_RES_KEYWORDS:
        msg = f"Output res: {sanitized_output_res} not one of {VALID_OUTPUT_RES_KEYWORDS}."
        raise ValueError(msg)


def apply_dn_mask(
    pvis: np.ndarray,
    vis_disp: np.ndarray,
    time_info_dt: dt.datetime,
    rlats: np.ndarray,
    rlons: np.ndarray,
) -> np.ndarray:
    """Creates combined day/night image by applying a day/night mask.

    Generates a day/night mask, and combines the adjusted Vis array and the
    ProxyVis array into a combined GeoProxyVis image.

    The pvis and vis_disp arguments could be of either 0.5 km or 2 km resolution,
    but BOTH arrays MUST have the same resolution.

    Args:
        pvis (np.ndarray): Full disk nighttime ProxyVis.
        vis (np.ndarray): Full disk input GEO Vis data.
        time_info_dt (dt.datetime): Time of the MIDDLE of the scan. This is NOT
            the GEO satellite data timestamp. The GEO timestamp is at the start
            of the scan. This is the time of the middle of the scan to calculate SZA most
            relevant for the full disk. For example, for GOES-16/18, and Himawari8/9 
            10-min full disk scan this is teh timestamp in the filename +5 minutes.
        rlats (np.ndarray): Array of lats matching resolution for pvis and vis_disp
        rlons (np.ndarray): Array of lons matching resolution for pvis and vis_disp

    Returns:
        pvis_combined (np.ndarray): Combined day/night data with adjusted Vis.
    """
    # start with 2 arrays both regridded to either 0.5 km or 2 km
    # nightime pvis, and daytime vis
    # apply dn_mask and make a single merged dn image

    # calculate SZA and day and night masks
    rsun_zen = solar_zenith_angle.sza_pysolar(time_info_dt, rlons, rlats)
    day_mask = solar_zenith_angle.sza_to_mask(
        np.float16(rsun_zen), mask_night=True, sza_thresh=SZA_THRESHOLD
    )
    night_mask = solar_zenith_angle.sza_to_mask(
        np.float16(rsun_zen), mask_night=False, sza_thresh=SZA_THRESHOLD
    )

    # make combined array
    pvis_combined = np.empty_like(pvis)

    night_mask_bool = np.isfinite(night_mask)
    day_mask_bool = np.isfinite(day_mask)

    pvis_combined[night_mask_bool] = VIS_SCALING_FACTOR * pvis[night_mask_bool]
    pvis_combined[day_mask_bool] = vis_disp[day_mask_bool]

    pvis_combined[~np.isfinite(pvis_combined)] = np.nanmax(pvis_combined)

    return pvis_combined


def calculate_pvis(
    satellite: str,
    data: DataDict,
    data_to_args: Dict[str, str],
    proxy_vis_alg_name: str,
    use_saved_params: bool,
) -> Tuple[np.ndarray, float, float]:
    """Performs user requested ProxyVis generation function.

    Args:
        satellite (str): Name of the satellite ("goes16", "goes17", "goes18",
            "himawari8", "himawari9", "meteosat-9", "meteosat-11")
        data (DataDict): Dictionary of IR and Vis data.
        data_to_args (Dict[str, str]): Dict mapping DataDict entries
            to pvis function args. Example: {'B07': 'c07', 'B11': 'c11'} 
        proxy_vis_alg_name (str): User requested ProxyVis function name:
            ("nighttime_pvis_main_two_eq", "nighttime_pvis_main_one_eq",
            "nighttime_pvis_simple_two_eq", "nighttime_pvis_simple_one_eq")
        use_saved_params (bool): Whether to use static/dynamic normalization
            for ProxyVis. Dynamic normalization works best when estimated from the full
            disk GEO data. For sub-sector ProxyVis images it is recommended to use saved
            normalization parameters.

    Returns:
        pvis_2km (np.ndarray): ProxyVis data at IR (2km) resolution.
        pvismin (float): The minimum value used to normalize the ProxyVis data.
            if use_saved_params was True, then this value was looked up based
            on the satellite.  Otherwise, this is calculated based on the data
            in data_dict.
        pvismax (float): The maximum value used to normalize the ProxyVis data.
            if use_saved_params was True, then this value was looked up based
            on the satellite.  Otherwise, this is calculated based on the data
            in data_dict.
    """
    apv = PVIS_FUNC_LOOKUP[proxy_vis_alg_name]

    channel_args = _create_channel_call_args(data, data_to_args, "bt_temp")

    # estimate ProxyVis
    pvis_2km, _pvis_regr, pvismin, pvismax = apv(
        satellite, **channel_args, use_saved_params=use_saved_params
    )

    # Only main pvis array is returned since other shoud not be needed for real-time processing
    return pvis_2km, pvismin, pvismax


def calculate_vis(
    satellite: str,
    data: DataDict,
    dvis_data_to_args: Dict[str, str],
    dvis_alg_name: str,
    rlons: np.ndarray,
    rlats: np.ndarray,
    time_info_dt: dt.datetime,
) -> Tuple[np.ndarray, float, float]:
    """Performs user requested Vis generation function.

    Args:
        satellite (str): Name of the satellite ("goes16", "goes17", "goes18",
            "himawari8", "himawari9", "meteosat-9", "meteosat-11")
        data (DataDict): Dictionary of IR and Vis data.
        dvis_data_to_args (Dict[str, str]): Dict mapping data entries
            to Vis function args. EX: {'B07': 'c07', 'B11': 'c11'}
        dvis_alg_name (str): Currently must be "vis_disp_sza"
        rlons (np.ndarray): Array of lons matching resolution for Vis data.
        rlats (np.ndarray): Array of lats matching resolution for Vis data.
        time_info_dt (dt.datetime): Time of the MIDDLE of the scan. This is NOT
            the GEO satellite data timestamp. The GEO timestamp is at the start
            of the scan. This is the time of the middle of the scan to calculate SZA most
            relevant for the full disk. For example, for GOES-16/18, and Himawari8/9 
            10-min full disk scan this is teh timestamp in the filename +5 minutes.

    Returns:
        vis_05km (np.ndarray): Adjusted Vis data.
        vismin (float): The minimum value used to normalize the Vis data.
        vismax (float): The maximum value used to normalize the Vis data.
    """
    apv = VIS_FUNC_LOOKUP[dvis_alg_name]

    channel_args = _create_channel_call_args(data, dvis_data_to_args, "radiances")

    # Generate adjusted Vis data
    vis_05km, vismin, vismax = apv(
        satellite, rlons, rlats, time_info_dt, **channel_args
    )

    return vis_05km, vismin, vismax


def _create_channel_call_args(
    data: DataDict,
    data_to_args: Dict[str, str],
    data_type: str,
) -> Dict[str, np.ndarray]:
    """Create Vis and ProxyVis func calling args using map from data dict names
    to call args. This is needed so that data from other geostationary satellites
    can be used with the subroutines that use ABI Vis and IR channels

    Creates a dictionary containing the channel calling args for ProxyVis and
    Vis functions. Uses a map to determine what entries in the data dict corespond
    to each channel argument in the user requested ProxyVis function.  For example,
    the Himawari AHI data dict contains "B02", "B07", "B11", "B13", and "B15"
    entries.  This needs to map to "c02" for Vis, and "c07", "c11", "c13", and
    "c15" for ProxyVis calling arguments.  The data_to_args map to accomplish this
    would be: {"B03":"c02} and {"B07":"c07, "B11":"c11", "B13":"c13", "B15":"c15"},
    for Vis and ProxyVis  correspondingly. The returned dictionary's keys are the
    function argument names and the values are the arrays corresponding to that
    argument.  The ** calling syntax can be used with the returned dictionary to
    pass the appropriate arrays to the user requested ProxyVis and Vis functions.

    Args:
        data (DataDict): Dictionary containing input GEO channel data.
        data_to_args (Dict[str, str]): Dictionary mapping contents of data
            dictionary to ProxyVis or Vis functions arguments.
        data_type: a string, 'bt_temp' for ProxyVis or 'readiances' for Vis

    Returns:
        args (Dict[str, np.ndarray]): Dictionary of data to be used as the
            channel arguments when calling a ProxyVis function.
    """
    args = {arg: data[data_type][chan] for chan, arg in data_to_args.items()}
    return args
