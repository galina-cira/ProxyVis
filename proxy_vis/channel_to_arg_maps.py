"""
    Convenience module providing maps from satellite channels to function input.

    The following module provides constants that can be used as the data_to_args_map
    and vis_channel arguments of the combine_dn_pvis.get_all_vis_pvis function.
    These map the contents of the data_dict to the arguments of the desired ProxyVis
    function. The vis channel that should be used is also provided. Arguments
    are provided for every combination of ProxyVis function and supported
    satellites.  The names used here are the names provided by Satpy.  If you are
    using a different reader, then you will need to make your own mapping with the
    names provided by your reader.

    The contstants with "MAIN" in their name will map channels to the arguments of
    the ProxyVis "main" functions.  The constants with "SIMPLE" in their name will
    map channels to the arguments of the ProxyVis "simple" functions.

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
    ##########################################################################
"""

# ABI GOES 16/17/18
ABI_MAIN = {"C07": "c07", "C11": "c11", "C13": "c13", "C15": "c15"}
ABI_SIMPLE = {"C07": "c07"}
ABI_VIS = {"C02": "c02"}

# AHI Himawari 8/9
AHI_MAIN = {"B07": "c07", "B11": "c11", "B13": "c13", "B15": "c15"}
AHI_SIMPLE = {"B07": "c07"}
AHI_VIS = {"B03": "c02"}

# SEVIRI Meteosat 2nd generation (MTG)
SEVIRI_MAIN = {
    "IR_039": "c07",
    "IR_087": "c11",
    "IR_108": "c13",
    "IR_120": "c15",
}
SEVIRI_SIMPLE = {"IR_039": "c07"}
SEVIRI_VIS = {"VIS006": "c02"}
