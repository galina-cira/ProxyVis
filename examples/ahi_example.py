"""
    ##########################################################################
    This is an example driver script that demonstrates how to use the high level
    ProxyVis API to process Himawari data.  Other satellites would be very
    similar.  This example assumes you have the proxy_vis, satpy, pyresample,
    and matplotlib libraries installed in your environment. To run, place B03,
    B07, B11, B13, and B15 Himawari-9 data for a single time in the
    examples/input/ directory.
    To run the code, you can use the command:
    python examples/ahi_example.py

    This will produce a plot of ProxyVis data in the examples/output/ directory.
    
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

import glob

import matplotlib.pyplot as plt
import numpy as np
import satpy

from proxy_vis import channel_to_arg_maps, combine_dn_pvis

# Use 32bit floats for the lons/lats generated from pyresample AreaDefinition
# objects.
LON_LAT_DTYPE = "f4"

# Satpy reads visible data in a range of 0 - 130, but the ProxyVis code assumes
# visible data is in a range of 0 - 1.3, so we need to scale this data.
SATPY_VIS_TO_PVIS_SCALE = 100.0

# Duration of each ahi scan in minutes.
AHI_SCAN_DURATION_MINUTES = 10

# Output ProxyVis composite provides data in range 0 - 1.3, so this is used to
# explicitly set the vmin/vmax parameters when plotting.
PLOT_VMIN = 0.0
PLOT_VMAX = 1.3

def main():
    # Read data using satpy
    filenames = glob.glob("examples/input/*.bz2")
    scene = satpy.Scene(filenames=filenames, reader="ahi_hsd")
    scene.load(["B03", "B07", "B11", "B13", "B15"])

    # Look up the start time from one of the channel's metadata
    start_time = scene["B07"].attrs["start_time"]

    # Convert the satpy Scene to a dictionary of numpy arrays
    # The dictionary uses a sub-dictionary to divide the channels between IR 
    # channels that use brightness temperatues ("bt_temp"), and visible channels
    # that use radiances("radiances").
    data_dict = {
        "bt_temp": {
            "B07": scene["B07"].data.compute(),
            "B11": scene["B11"].data.compute(),
            "B13": scene["B13"].data.compute(),
            "B15": scene["B15"].data.compute(),
        },
        "radiances": {"B03": scene["B03"].data.compute() / SATPY_VIS_TO_PVIS_SCALE},
    }

    # Compute the 2km and 0.5km lons/lats from the pyresample AreaDefinition
    # included in the satpy scene metadata.
    lons_2km, lats_2km = scene["B07"].attrs["area"].get_lonlats(dtype=LON_LAT_DTYPE)
    lons_0p5km, lats_0p5km = scene["B03"].attrs["area"].get_lonlats(dtype=LON_LAT_DTYPE)

    # We're just going to generate a 2km ProxyVis image.
    # _pvis_combined_05km data will be None.
    _pvis_combined_05km, pvis_combined_2km, _pvismin, _pvismax = (
        combine_dn_pvis.get_all_vis_pvis(
            "himawari9", 
            start_time,
            data_dict,
            channel_to_arg_maps.AHI_MAIN,
            channel_to_arg_maps.AHI_VIS,
            lons_2km,
            lats_2km,
            lons_0p5km,
            lats_0p5km,
            AHI_SCAN_DURATION_MINUTES,
            "nighttime_pvis_main_two_eq",
            "vis_disp_sza",
            True,
            "2.0km", # Just produce 2km data.
        )
    )

    # Set pixels off the disk to NaN by finding NaNs and INFs in longitude
    # array. Appears to work for all the satellites tested used with ProxyVis.
    # An alternative would be to calculate the max valid distance or the max
    # valid satellite zenith angle.
    off_disk = ~np.isfinite(lons_2km)
    pvis_combined_2km[off_disk] = np.nan

    output_filename = f"examples/output/ahi_2km_{start_time:%Y-%m-%dT%H}.png"
    plot(pvis_combined_2km, output_filename)


def plot(data: np.ndarray, filename: str) -> None:
    fig = plt.figure(figsize=(10, 10))
    fig.tight_layout()
    ax = fig.add_subplot(1, 1, 1)
    ax.axis("off")
    ax.imshow(data, vmin=PLOT_VMIN, vmax=PLOT_VMAX, cmap="Greys_r")
    fig.savefig(filename, dpi=550, bbox_inches="tight", pad_inches=0, transparent=True)

    fig.clf()


if __name__ == "__main__":
    main()
