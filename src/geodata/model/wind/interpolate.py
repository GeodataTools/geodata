# Copyright 2024 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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
import scipy.interpolate as sinterp
import xarray as xr

from ...logging import logger
from ._base import WindBaseModel

# See https://confluence.ecmwf.int/display/UDOC/L137+model+level+definitions
LEVEL_TO_HEIGHT = {
    131: 169.5,
    132: 136.62,
    133: 106.54,
    134: 79.04,
    135: 53.92,
    136: 30.96,
    137: 10.0,
}


class WindInterpolationModel(WindBaseModel):
    """Wind speed estimation based on the an interpolation model.
    Example:
        >>> from geodata import Dataset
        >>> from geodata.model.wind import WindInterpolationModel
        >>> dataset = Dataset(module="era5", weather_data_config="wind_3d_hourly", years=slice(2010, 2010), months=slice(1,2))
        >>> model = WindExtrapolationModel(dataset)
        >>> model.prepare()
        >>> model.estimate(height=12, xs=slice(1, 2), ys=slice(1, 2), years=slice(2010, 2010), months=slice(1, 2))
    """

    SUPPORTED_WEATHER_DATA_CONFIGS = {"wind_3d_hourly"}

    def _prepare_fn(
        self,
        ds: xr.Dataset,
        half_precision: bool = True,
    ) -> xr.Dataset:
        """Compute wind speed using the ERA5 3D dataset.

        Args:
            ds (xr.Dataset): ERA5 3D dataset.
            half_precision (bool, optional): Use float32 precision to store coefficients and residuals. Defaults to True.

        Returns:
            xr.Dataset: Dataset with wind speed.
        """

        assert (
            "model_level" in ds.coords
        ), "Dataset does not contain model levels. Please double-check the dataset."

        ds = ds.transpose("model_level", ...).sortby(
            "model_level", ascending=False
        )  # Model levels are sorted in descending order

        variables = [var for var in ds if var.replace("u", "v") in ds and var == "u"]
        heights = np.array(
            [LEVEL_TO_HEIGHT[int(level)] for level in ds["model_level"].values]
        )

        logger.debug("Selected variables: %s", variables)
        logger.debug("Shape of heights: %s", heights.shape)

        speeds = (ds["u"] ** 2 + ds["v"] ** 2) ** 0.5
        orig_shape = speeds.shape
        speeds = speeds.values.reshape(len(heights), -1)

        spline_params: sinterp.BSpline = sinterp.make_interp_spline(
            heights, speeds, k=3
        )
        t, c = spline_params.t, spline_params.c

        if half_precision:
            c = c.astype(np.float32)
        ds["c"] = ds["u"].copy(data=c.reshape(orig_shape))
        ds.attrs["t"] = t

        return ds[["c"]]

    def _estimate_dataset(
        self,
        height: int,
        years: slice,
        months: Optional[slice] = None,
        xs: Optional[slice] = None,
        ys: Optional[slice] = None,
        use_real_data: Optional[bool] = False,
    ) -> xr.Dataset:
        params = xr.open_mfdataset(self.files).transpose("model_level", ...)
        spline_params = sinterp.BSpline(
            params.attrs.get("t"), params.get("c").values, k=3
        )

        params = params.drop_dims("model_level")
        speeds = xr.DataArray(
            spline_params(height), dims=params.dims, coords=params.coords
        )

        return speeds

    def _estimate_cutout(
        self,
        height: int,
        years: slice,
        months: Optional[slice] = None,
        xs: Optional[slice] = None,
        ys: Optional[slice] = None,
        use_real_data: Optional[bool] = False,
    ) -> xr.Dataset:
        pass
