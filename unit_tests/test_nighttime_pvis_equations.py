""" unit tests for night time ProxyVis calculations 

Input and output arrays were generated with original code to verify that the new code produces the same results

"""
import pytest

import numpy as np

from proxy_vis.nighttime_pvis_main_one_eq import calculate_pvis_main_one_eq
from proxy_vis.nighttime_pvis_main_two_eq import calculate_pvis_main_two_eq
from proxy_vis.nighttime_pvis_simple_one_eq import calculate_pvis_simple_one_eq
from proxy_vis.nighttime_pvis_simple_two_eq import calculate_pvis_simple_two_eq

# pylint: disable=line-too-long
INPUT_DATA = {
        'c15': np.array([283.96, 211.83, 224.61, 278.67, 286.88, 259.49, 274.67, 276.41, 290.22, 286.98, 288.95, 284.67, 274.56, 253.62, 249.77, 203.43]),
        'c13': np.array([284.28, 216.89, 228.09, 281.85, 291.18, 267.85, 276.6, 279.14, 294.63, 291.36, 291.96, 288.85, 279.86, 256.3, 248.7, 205.65]),
        'c11': np.array([282.62, 217.18, 228.19, 279.14, 288.27, 267.51, 274.48, 275.72, 291.89, 288.36, 288.21, 283.63, 273.14, 252.62, 246.32, 204.79]),
        'c07': np.array([283.59, 232.19, 238.29, 281.43, 293.73, 281.09, 276.04, 280.14, 297.92, 295.26, 293.04, 290.02, 281.76, 264.82, 251.48, 197.31]),
        }

main_eq_one_g16 = np.array([0.550,0.791,0.783,0.510,0.304,0.426,0.588,0.508,0.220,0.267,0.333,0.359,0.450,0.622,0.751,0.920])
pvisible_simpleA_g16 = np.array([0.477,0.849,0.822,0.503,0.333,0.507,0.561,0.517,0.256,0.307,0.345,0.391,0.499,0.661,0.753,0.952])
pvisible_simple_g16 = np.array([0.507,0.836,0.810,0.538,0.326,0.543,0.610,0.556,0.224,0.291,0.340,0.400,0.534,0.650,0.741,0.938])
main_eq_two_g16 = np.array([0.571, 0.758, 0.756, 0.547, 0.308, 0.471, 0.632, 0.545, 0.204, 0.264, 0.335, 0.374, 0.489, 0.613, 0.741, 0.897])

main_eq_one_g17 = np.array([0.524, 0.753, 0.745, 0.486, 0.29, 0.406, 0.56, 0.484, 0.209, 0.254, 0.317, 0.342, 0.429, 0.593, 0.715, 0.882])
main_eq_two_g17 = np.array([0.543, 0.721, 0.719, 0.52, 0.293, 0.449, 0.601, 0.519, 0.194, 0.251, 0.319, 0.356, 0.465, 0.584, 0.706, 0.854])
pvisible_simple_g17 = np.array([0.483, 0.797, 0.772, 0.513, 0.31, 0.518, 0.581, 0.531, 0.213, 0.278, 0.324, 0.381, 0.509, 0.62, 0.706, 0.893])
pvisible_simpleA_g17 = np.array([0.455, 0.808, 0.783, 0.479, 0.317, 0.483, 0.534, 0.493, 0.244, 0.292, 0.328, 0.373, 0.475, 0.629, 0.717, 0.906])

main_eq_one_him = np.array([0.546, 0.784, 0.776, 0.506, 0.302, 0.423, 0.583, 0.505, 0.218, 0.265, 0.33, 0.356, 0.447, 0.617, 0.745, 0.919])
main_eq_two_him = np.array([0.566, 0.751, 0.749, 0.542, 0.305, 0.467, 0.627, 0.54, 0.202, 0.261, 0.333, 0.371, 0.484, 0.608, 0.735, 0.889])
pvisible_simple_him = np.array([0.504, 0.83, 0.804, 0.535, 0.323, 0.54, 0.606, 0.553, 0.222, 0.289, 0.338, 0.397, 0.53, 0.646, 0.736, 0.931])
pvisible_simpleA_him = np.array([0.474, 0.842, 0.816, 0.499, 0.331, 0.503, 0.556, 0.513, 0.254, 0.304, 0.342, 0.388, 0.495, 0.656, 0.746, 0.944])
#pylint: enable=line-too-long


tests = [(calculate_pvis_main_one_eq, main_eq_one_g16, 'goes16', 0.78),
         (calculate_pvis_main_two_eq, main_eq_two_g16, 'goes16', 0.78),
         (calculate_pvis_main_one_eq, main_eq_one_g17, 'goes17', 0.84),
         (calculate_pvis_main_two_eq, main_eq_two_g17, 'goes17', 0.84),
         (calculate_pvis_main_one_eq, main_eq_one_him, 'himawari9', 0.79),
         (calculate_pvis_main_two_eq, main_eq_two_him, 'himawari9', 0.79),
         ]
@pytest.mark.parametrize("func, return_array, satellite, exp_prmax", tests)
def test_calculate_prvis_main(func, return_array, satellite, exp_prmax):
    """ Test case multi-channel equations """
    proxy_vis, _, prvismin, prvismax = func(satellite, **INPUT_DATA, use_saved_params=True)
    assert np.isclose(prvismin, 0.0, rtol=1e-5, atol=1e-5)
    assert np.isclose(prvismax, exp_prmax, rtol=1e-5, atol=1e-5)
    assert pytest.approx(return_array, rel=1e-2, abs=1e-2) ==  proxy_vis


tests = [
         (calculate_pvis_simple_one_eq, pvisible_simpleA_g16, 'goes16', 0.78),
         (calculate_pvis_simple_two_eq, pvisible_simple_g16, 'goes16', 0.78),
         (calculate_pvis_simple_one_eq, pvisible_simpleA_g17, 'goes17', 0.84),
         (calculate_pvis_simple_two_eq, pvisible_simple_g17, 'goes17', 0.84),
         (calculate_pvis_simple_one_eq, pvisible_simpleA_him, 'himawari9', 0.79),
         (calculate_pvis_simple_two_eq, pvisible_simple_him, 'himawari9', 0.79),
         ]
@pytest.mark.parametrize("func, return_array, satellite, exp_prmax", tests)
def test_calculate_prvis_simple(func, return_array, satellite, exp_prmax):
    """ Test casesimple, single channel equations """
    proxy_vis, _, prvismin, prvismax = func(satellite, INPUT_DATA['c07'], use_saved_params=True)
    assert np.isclose(prvismin, 0.0, rtol=1e-5, atol=1e-5)
    assert np.isclose(prvismax, exp_prmax, rtol=1e-5, atol=1e-5)
    assert pytest.approx(return_array, rel=1e-2, abs=1e-2) ==  proxy_vis
