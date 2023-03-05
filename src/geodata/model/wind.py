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

import numpy as np
import xarray as xr
from numba import njit, prange

logger = logging.getLogger(__name__)


HEIGHTS = {"u50m": 50, "u10m": 10, "u2m": 2}


@njit(parallel=True)
def _compute_wind_speed(
    heights: np.ndarray, speeds: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Compute wind speed from heights and speeds.

    Args:
        heights (np.ndarray): Array of heights.
        speeds (np.ndarray): Array of speeds.

    Returns:
        tuple[np.ndarray, np.ndarray]: Tuple of coefficients (alpha, beta) and residuals.
    """
    # pylint: disable=not-an-iterable
    coeffs = np.empty(shape=heights.shape[:-1] + (2,))
    residuals = np.empty_like(speeds)
    for time in prange(heights.shape[0]):
        for lat in prange(heights.shape[1]):
            for lon in prange(heights.shape[2]):
                a = heights[time, lat, lon]
                a = np.stack((a, np.ones_like(a)), axis=-1)
                b = np.abs(speeds[time, lat, lon])

                coeff, residual, _, _ = np.linalg.lstsq(a, b)

                coeffs[time, lat, lon] = coeff
                residuals[time, lat, lon] = residual
    # pylint: enable=not-an-iterable

    return coeffs, residuals


def wind_speed_merra2(ds: xr.Dataset, compute_lml: bool = True) -> xr.Dataset:
    """Compute wind speed from for MERRA2 dataset.

    Args:
        ds (xr.Dataset): MERRA2 dataset.
        compute_lml (bool, optional): Include speed at LML in calculation. Defaults to True.

    Returns:
        xr.Dataset: Dataset with wind speed.
    """

    disph = ds["disph"].values

    variables = [f for f in HEIGHTS if f in ds and f.replace("u", "v") in ds]
    heights = np.array([HEIGHTS[f] for f in variables]) - disph[..., np.newaxis]

    if compute_lml:
        hlml = ds["hlml"].values
        variables.append("ulml")
        heights = np.concatenate((heights, (hlml - disph)[..., np.newaxis]), axis=-1)
    log_heights = np.log(heights, out=-np.ones_like(heights), where=heights > 0)

    speeds = []
    for var in variables:
        speeds.append(((ds[var] ** 2 + ds[var.replace("u", "v")] ** 2) ** 0.5).values)
    speeds = np.stack(speeds, axis=-1).astype(heights.dtype)

    coeffs, residuals = _compute_wind_speed(log_heights, speeds)
    ds["coeffs"] = (("time", "lat", "lon", "coeff"), coeffs)
    ds["residuals"] = (("time", "lat", "lon", "residual"), residuals)

    return ds
