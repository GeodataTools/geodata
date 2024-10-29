# Copyright 2023 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from typing import Optional

import numpy as np
import pandas as pd
import xarray as xr

from ...logging import logger
from ._base import HEIGHTS, WindBaseModel

try:
    from numba import njit, prange
except ImportError:
    logger.warning("Numba not installed. Using pure Python implementation.")
    prange = range

    from ...utils import dummy_njit as njit


@njit(parallel=True)
def _compute_wind_speed_ext(
    heights: np.ndarray, speeds: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Compute wind speed from heights and speeds.

    Args:
        heights (np.ndarray): Array of heights.
        speeds (np.ndarray): Array of speeds.

    Returns:
        tuple[np.ndarray, np.ndarray]: Tuple of coefficients (alpha, beta) and residuals.
    """

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

    return coeffs, residuals


class WindExtrapolationModel(WindBaseModel):
    """Wind speed estimation based on the an extrapolation model.

    Model Details:
        The wind speed is estimated using a linear regression of the logarithm of the wind speed
        against the logarithm of the height. The coefficients of the regression are stored in the
        dataset as a 2D array of shape (2,). The first coefficient is the slope of the regression
        line, and the second coefficient is the intercept. The wind speed is then estimated as
        :math:`\\hat{v} = \\alpha \\cdot \\log(z) + \\beta`, where :math:`\\alpha` and :math:`\\beta`
        are the coefficients, and :math:`z` is the height.

    Residuals:
        The residuals of the regression are stored in the dataset as a 2D array of shape (n,).
        The residuals are the difference between the estimated wind speed and the actual wind speed.
        The residuals are stored in the same order as the heights used in the regression.

    Example:
        >>> from geodata import Dataset
        >>> from geodata.model.wind import WindExtrapolationModel
        >>> dataset = Dataset(module="merra2", weather_data_config="slv_flux_hourly", years=slice(2010, 2010), months=slice(1,2))
        >>> model = WindExtrapolationModel(dataset)
        >>> model.prepare()
        >>> model.estimate(height=12, xs=slice(1, 2), ys=slice(1, 2), years=slice(2010, 2010), months=slice(1, 2))
    """

    SUPPORTED_WEATHER_DATA_CONFIGS = {"slv_flux_hourly"}

    def _prepare_fn(
        self,
        ds: xr.Dataset,
        compute_lml: bool = True,
        half_precision: bool = True,
        compute_residuals: bool = True,
    ) -> xr.Dataset:
        """Compute wind speed from for MERRA2 dataset.

        Args:
            ds (xr.Dataset): MERRA2 dataset.
            compute_lml (bool, optional): Include speed at LML in calculation. Defaults to True.
            half_precision (bool, optional): Use float32 precision to store coefficients and residuals. Defaults to True.
            compute_residuals (bool, optional): Compute and store residuals against existing data points. Defaults to True.

        Returns:
            xr.Dataset: Dataset with wind speed.
        """

        disph = ds["disph"].values

        variables = [f for f in HEIGHTS if f in ds and f.replace("u", "v") in ds]
        heights = np.array([HEIGHTS[f] for f in variables]) - disph[..., np.newaxis]

        logger.debug("Selected variables: %s", variables)
        logger.debug("Shape of heights: %s", heights.shape)
        # Basic sanity check
        if 0 in heights.shape:
            raise ValueError(
                "Dataset does not contain any other useable heights other than lml"
            )

        if compute_lml:
            hlml = ds["hlml"].values
            variables.append("ulml")
            heights = np.concatenate(
                (heights, (hlml - disph)[..., np.newaxis]), axis=-1
            )
        log_heights = np.log(heights, out=-np.ones_like(heights), where=heights > 0)

        speeds = []
        for var in variables:
            speeds.append(
                ((ds[var] ** 2 + ds[var.replace("u", "v")] ** 2) ** 0.5).values
            )
        speeds = np.stack(speeds, axis=-1).astype(heights.dtype)

        coeffs, residuals = _compute_wind_speed_ext(log_heights, speeds)

        if half_precision:
            coeffs = coeffs.astype("float32")
            residuals = residuals.astype("float32")
        ds = ds.assign_coords(coeff=["alpha", "beta"])
        ds["coeffs"] = (("time", "lat", "lon", "coeff"), coeffs)

        if compute_residuals:
            ds = ds.assign_coords(residual=variables)
            ds["residuals"] = (("time", "lat", "lon", "residual"), residuals)

        return ds[["coeffs", "residuals"]] if compute_residuals else ds[["coeffs"]]

    def _estimate_dataset(
        self,
        height: int,
        years: slice,
        months: Optional[slice] = None,
        xs: Optional[slice] = None,
        ys: Optional[slice] = None,
        use_real_data: Optional[bool] = False,
    ) -> xr.DataArray:
        assert height > 0, "Height must be greater than 0."

        if months is None:
            months = slice(1, 12)

        start_time = pd.Timestamp(year=years.start, month=months.start, day=1)
        end_time = pd.Timestamp(
            year=years.stop, month=months.stop, day=31, hour=23, minute=59, second=59
        )

        ds = xr.open_mfdataset(self.files)

        if xs is None:
            xs = ds.coords["lon"]
        if ys is None:
            ys = ds.coords["lat"]

        ds = ds.sel(lon=xs, lat=ys, time=slice(start_time, end_time))

        if height in HEIGHTS.values() and use_real_data:
            logger.info("Using real data for estimation at height %d", height)
            return (ds[f"u{height}m"] ** 2 + ds[f"v{height}m"] ** 2) ** 0.5

        alpha = ds["coeffs"][..., 0]
        beta = ds["coeffs"][..., 1]

        result = alpha * np.log((height - ds["disph"]) / np.exp(-beta / alpha))
        return result.drop_vars("coeff")  # remove unnecessary coordinate

    def _estimate_cutout(
        self,
        height: int,
        years: slice,
        months: Optional[slice] = None,
        xs: Optional[slice] = None,
        ys: Optional[slice] = None,
        use_real_data: Optional[bool] = False,
    ) -> xr.DataArray:
        assert height > 0, "Height must be greater than 0."

        if months is None:
            months = slice(1, 12)

        start_time = pd.Timestamp(year=years.start, month=months.start, day=1)
        end_time = pd.Timestamp(
            year=years.stop, month=months.stop, day=31, hour=23, minute=59, second=59
        )

        ds = xr.open_mfdataset(self.files)

        if xs is None:
            xs = ds.coords["x"]
        if ys is None:
            ys = ds.coords["y"]

        ds = ds.sel(x=xs, y=ys, time=slice(start_time, end_time))

        if height in HEIGHTS.values() and use_real_data:
            logger.info("Using real data for estimation at height %d", height)
            return (ds[f"u{height}m"] ** 2 + ds[f"v{height}m"] ** 2) ** 0.5

        alpha = ds["coeffs"][..., 0]
        beta = ds["coeffs"][..., 1]

        result = alpha * np.log((height - ds["disph"]) / np.exp(-beta / alpha))
        return result.drop_vars("coeff")  # remove unnecessary coordinate
