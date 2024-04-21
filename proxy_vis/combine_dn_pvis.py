"""
    This is the main module to esimate both daytime Vis and nighttime ProxyVis
    parts of the GeoProxyVis imagery your code shoud call the soubroutine
    get_all_vis_pvis

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
from typing import Dict

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

DataDict = Dict[str, Dict[str, np.ndarray]]
"""Dictionary whose keys are channel names and values are dictionaries of data.

Dictionary of data read from a file.  The keys should be channel names and the
values should be a sub-dictionary of calibrated data.  The keys for the 
sub-dictionary should be the type of calibration, one of "bt_temp", "radiances".
"""

# Available ProxyVis functions
PVIS_FUNC_LOOKUP = {
    "nighttime_pvis_main_two_eq": nighttime_pvis_main_two_eq.calculate_pvis_main_two_eq,
    "nighttime_pvis_main_one_eq": nighttime_pvis_main_one_eq.calculate_pvis_main_one_eq,
    "nighttime_pvis_simple_two_eq": nighttime_pvis_simple_two_eq.calculate_pvis_simple_two_eq,
    "nighttime_pvis_simple_one_eq": nighttime_pvis_simple_one_eq.calculate_pvis_simple_one_eq,
}

# Available Visible functions
VIS_FUNC_LOOKUP = {"vis_disp_sza": daytime_pvis.vis_disp_sza}

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
):
    """
    Get combined day/time full disk GeoProxyVis imagery

    inputs:
        satellite       :   string, satellite ID (goes-16, -goes-17, himawari)
        channels        :   list; list of input channels
        ch_letter       :   satellite specific letter for channel
                                ( "C" for GOES-16/17, "B" for Himawari )

        time_dt         :   datetime object
                            time of the HW full disk image. The time *MUST* match
                                    exactly the timestamp in HW filename.
        data_dict       :  dictionary: keys - HW channels ('B03', etc),
                            values - corresponding data arrays. This must be verion
                                with all bad values masked
        geo_lons_05km   : numpy array
                                array of lons for HW naviagtion for 0.5 km version, numpy array
        geo_lats_05km   :  numpy array
                               array of lats for HW naviagtion for 0.5 km version, numpy array
        geo_lons_2km    :   numpy array
                                array of lons for HW naviagtion for 2 km version, numpy array
        geo_lats_2km    :   numpy array
                                array of lats for HW naviagtion for 2 km version, numpy array
        pvis_alg   :  function to use to calculate proxy_vis
        dvis_alg   :  function to use to calculate daytime vis

        use_saved_params    : boolean;
                            : True  : use pre-saved min/max for ProxyVis
                            : False : estimates min/max for ProxyVis from data at run time
                                Note: min/max MUST be estimated from the full disk. If processing
                                is fo a subsection of disk, pre-saved parameters MUST be used to get
                                correct output

    outputs:

        pvis_combined_05km  :   numpy array
                                    full disk combined day/night ProxyVis at 0.5 km resolution; numpy array
        pvismin            :   float
                                    min values used for proxy_vis
        pvismax            :   float
                                    max value used for proxy_vis


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
):
    """
    Apply day/night mask to Vis and Pvis arrays and create combined arrays
    Also creates adjusted Vis array

    inputs:
        pvis        : full disk nighttime ProxyVis; numpy array
        vis         : full disk non-adjusted original vis: numpy array

        NOTE: pvis and vis could be of either 0.5 km or 2 km resolution,
        however BOTH arrays MUTS have the same resolution

        time_info_dt    : datetime object
                        : time of the MIDDLE of the scan. This is NOT the HW timestamp.
                            The HW timestamp is at the start of the scan. This is
                            at the middle of the scan to calculate SZA most relevant for th efull disk.

        rlats       : numpy arrys of lats matching resolution for pvis and vis
        rlons       : numpy arrys of lons matching resolution for pvis and vis

    outputs:

        pvis_combined   : numpy array
                        : combinde day/night data with adjusted Vis

    """
    # start with 2 arrays both regridded to eother 0.5 km or 2 km
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
):
    """
    Prepare data to call the subroutine to calculate nighttime part of
    GeoProxyVis


    inputs:
        satellite       :   string, satellite ID (goes16, goes17, himawari)
        data:           :   dictionary
                        :   keys - HW channels ('B03','B07', etc), values - corresponding data arrays
                                returned is version with allvalues outside thresh_dist set to NaN
        data_to_args    :   Dictionary mapping the names of bands in the data
                            dictionary to the calling arguments of the desired
                            ProxyVis function.  EX: {'B07': 'c07', 'B11': 'c11'}
        proxy_vis_alg   :   function name to use to calculate proxy_vis, one of
                            nighttime_pvis_main_two_eq, nighttime_pvis_main_one_eq,
                            nighttime_pvis_simple_two_eq, nighttime_pvis_simple_one_eq
        use_saved_params    : boolean;
                            : True  : use pre-saved min/max for ProxyVis
                            : False : estimates min/max for ProxyVis from data at run time
                                Note: min/max MUST be estimated from the full disk. If processing
                                is fo a subsection of disk, pre-saved parameters MUST be used to get
                                correct output


    outputs:
        pvis_2km            : full disk nighttime ProxyVis at original 2km resolution; numpy array
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
):
    """
    Prepare data to call the subroutine to calculate daytime part of
    GeoProxyVis


    inputs:
        satellite       :   string, satellite ID (goes16, goes17, himawari)
        data:           :   dictionary
                        :   keys - HW channels ('B03','B07', etc), values - corresponding data arrays
                                returned is version with allvalues outside thresh_dist set to NaN
        data_to_args    :   Dictionary mapping the names of bands in the data
                            dictionary to the calling arguments of the desired
                            ProxyVis function.  EX: {'B07': 'c07', 'B11': 'c11'}
        proxy_vis_alg   :   function name to use to calculate proxy_vis, one of
                            nighttime_pvis_main_two_eq, nighttime_pvis_main_one_eq,
                            nighttime_pvis_simple_two_eq, nighttime_pvis_simple_one_eq
        use_saved_params    : boolean;
                            : True  : use pre-saved min/max for ProxyVis
                            : False : estimates min/max for ProxyVis from data at run time
                                Note: min/max MUST be estimated from the full disk. If processing
                                is fo a subsection of disk, pre-saved parameters MUST be used to get
                                correct output


    outputs:
        pvis_2km            : full disk nighttime ProxyVis at original 2km resolution; numpy array
    """
    apv = VIS_FUNC_LOOKUP[dvis_alg_name]

    channel_args = _create_channel_call_args(data, dvis_data_to_args, "radiances")

    # estimate ProxyVis
    vis_05km, vismin, vismax = apv(
        satellite, rlons, rlats, time_info_dt, **channel_args
    )

    # Only main pvis array is returned since other shoud not be needed for real-time processing
    return vis_05km, vismin, vismax


def _create_channel_call_args(
    data: DataDict,
    data_to_args: Dict[str, str],
    data_type: str,
) -> Dict[str, np.ndarray]:
    """Create Vis and ProxyVis func calling args using map from data dict names
    to call args.  This is needed so that data from other geostationary satellites
    can be used with th esubroutines that use ABI Vis and IR channels

    Creates a dictionary containing the channel calling args for a ProxyVis and Vis
    functions. Uses a map to determine what entries in the data dict corespond to
    each channel argument in the desired ProxyVis function.  For example, the Himawari
    AHI data dict contains i"B02", "B07", "B11", "B13", and "B15" entries.  This
    needs to map to "c02" for Vis, and "c07", "c11", "c13", and "c15" for ProxyVis calling arguments.  The
    data_to_args map to accomplish this would be: {"B03":"c02} and
    {"B07":"c07, "B11":"c11", "B13":"c13", "B15":"c15"}, fopr Vis and proxyVis  correspondingly.The returned
    dictionary's keys are the function argument names and the values are the
    arrays corresponding to that argument.  The ** calling syntax can be used
    with the returned dictionary to pass the appropriate arrays to the desired
    ProxyVis and Vis functions.

    Args:
        data (DataDict): Dictionary containing input GEO channel data.
        data_to_args (Dict[str, str]): Dictionary mapping contents of data
            dictionary to ProxyVis or Vis functions arguments.
        data_type: a string, 'bt_temp' for ProxyVis or 'readiances' for Vis

    Returns:
        Dict[str, np.ndarray]: Dictionary of data to be used as the channel
            arguments when calling a ProxyVis function.
    """
    args = {arg: data[data_type][chan] for chan, arg in data_to_args.items()}
    return args
