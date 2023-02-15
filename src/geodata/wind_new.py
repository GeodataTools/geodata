## Copyright 2023 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

import logging
import warnings

import numpy as np
import xarray as xr
from numba import njit, prange

import geodata

logger = logging.getLogger(__name__)


U_HEIGHTS = {"u50m": 50, "u10m": 10, "u2m": 2}
V_HEIGHTS = {"v50m": 50, "v10m": 10, "v2m": 2}


@njit(parallel=True)
def _compute_wind_speed(
    heights: np.ndarray, speeds: np.ndarray, return_residuals: bool = False
):
    """Compute wind speed from heights and speeds.

    Parameters
    ----------
    heights : np.ndarray
        Heights of wind speed measurements. It has Shape (time, lat, lon, n_heights).
    speeds : np.ndarray
        Wind speed measurements.

    Returns
    -------
    np.ndarray
        Wind speed.
    """

    # pylint: disable=not-an-iterable
    coeffs = np.empty(shape=heights.shape[:-1] + (2,))
    residuals = np.empty_like(speeds)
    for time in prange(heights.shape[0]):
        for lat in prange(heights.shape[1]):
            for lon in prange(heights.shape[2]):
                coeff, residual, _, _ = np.linalg.lstsq(
                    heights[time, lat, lon], speeds[time, lat, lon]
                )
                coeffs[time, lat, lon] = coeff
                residuals[time, lat, lon] = residual
    # pylint: enable=not-an-iterable

    return coeffs, residuals if return_residuals else coeffs


def wind_speed_merra2(
    ds: xr.Dataset, return_u: bool = True, return_v: bool = True
) -> xr.Dataset:
    """Compute wind speed from for MERRA2 dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing u and v components of wind and displacement height (disph).
    return_u : bool, optional
        Whether to compute u component of wind, by default True
    return_v : bool, optional
        Whether to compute v component of wind, by default True

    Returns
    -------
    xr.Dataset
        Dataset containing wind speed.
    """

    disph = ds["disph"].values

    if return_u:
        variables = [f for f in U_HEIGHTS if f in ds]
        heights = np.array([U_HEIGHTS[f] for f in variables]) - disph[..., np.newaxis]

        with warnings.catch_warnings():
            log_heights = np.concatenate(
                (np.ones_like(heights[..., 0]), np.log(heights)), axis=-1
            )

        speeds = []
        for var in variables:
            speeds.append(ds[var].values)
        speeds = np.stack(speeds, axis=-1)

        coeffs = _compute_wind_speed(log_heights, speeds)
        ds["u_coeffs"] = (("time", "lat", "lon"), coeffs)

    if return_v:
        variables = [f for f in V_HEIGHTS if f in ds]
        heights = np.array([V_HEIGHTS[f] for f in variables]) - disph[..., np.newaxis]

        with warnings.catch_warnings():
            log_heights = np.concatenate(
                (np.ones_like(heights[..., 0]), np.log(heights)), axis=-1
            )

        speeds = []
        for var in variables:
            speeds.append(ds[var].values)
        speeds = np.stack(speeds, axis=-1)

        coeffs = _compute_wind_speed(log_heights, speeds)
        ds["v_coeffs"] = (("time", "lat", "lon"), coeffs)

    return ds
