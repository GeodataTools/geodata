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

from typing import Hashable, Optional

import numpy as np
import scipy.interpolate as sinterp
import xarray as xr
from xarray.namedarray.pycompat import array_type as dask_array_type

from ...logging import logger
from ...utils import get_daterange
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


def _memoryview_safe(x: np.ndarray) -> np.ndarray:
    """Make array safe to run in a Cython memoryview-based kernel. These
    kernels typically break down with the error ``ValueError: buffer source
    array is read-only`` when running in dask distributed.

    Borrowed from https://github.com/crusaderky/xarray_extras/blob/main/xarray_extras/kernels/interpolate.py
    """
    if not x.flags.writeable:
        if not x.flags.owndata:
            x = x.copy(order="C")
        x.setflags(write=True)
    return x


def _make_interp_coeff(*args, **kwargs):
    """Dummy function to handle interpolation coefficients."""
    return sinterp.make_interp_spline(*args, **kwargs).c


def _splrep(a: xr.DataArray, dim: Hashable, k: int = 3) -> xr.Dataset:
    """Modified version of scipy.interpolate.splrep for xarray DataArray.
    Borrowed from https://github.com/crusaderky/xarray_extras/blob/main/xarray_extras/kernels/interpolate.py

    Args:
        a (xr.DataArray): Input data array.
        dim (Hashable): Dimension to interpolate along.
        k (int, optional): Degree of the spline. Defaults to 3.

    Returns:
        xr.Dataset: Dataset containing spline parameters.
    """

    # Make sure that dim is on axis 0
    a = a.transpose(dim, ...)
    x: np.ndarray = a.coords[dim].values

    if x.dtype.kind == "M":
        # Same treatment will be applied to x_new.
        # Allow x_new.dtype==M8[D] and x.dtype==M8[ns], or vice versa
        x = x.astype("M8[ns]").astype(float)

    t = sinterp._bsplines._not_a_knot(x, k=k)

    if isinstance(a.data, dask_array_type("dask")):
        from dask.array import map_blocks

        if len(a.data.chunks[0]) > 1:
            raise NotImplementedError(
                "Unsupported: multiple chunks on interpolation dim"
            )

        c = map_blocks(
            _make_interp_coeff,
            x,
            a.data,
            k=k,
            t=t,
            check_finite=False,
            dtype=float,
        )
    else:
        c = _make_interp_coeff(x, a.data, k=k, t=t, check_finite=False)

    return xr.Dataset(
        data_vars={
            "c": (a.dims, c),
        },
        coords=a.coords,
        attrs={
            "spline_dim": dim,
            "k": k,
            "t": t,
        },
    )


class WindInterpolationModel(WindBaseModel):
    """Wind speed estimation based on a spline interpolation of the wind speed at different heights.

    This model uses the ERA5 3D dataset to estimate wind speed at a given height using
    spline interpolation.

    Example:

    >>> from geodata import Dataset
    >>> from geodata.model.wind import WindInterpolationModel
    >>> dataset = Dataset(module="era5", weather_data_config="wind_3d_hourly", years=slice(2010, 2010), months=slice(1, 2))
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

        ds.coords["model_level"] = np.array(
            [LEVEL_TO_HEIGHT[int(level)] for level in ds["model_level"].values]
        )
        ds = (
            ds.rename({"model_level": "height"})
            .transpose("height", ...)
            .sortby("height")
        )

        logger.debug("Shape of heights: %s", ds["height"].shape)
        speeds = (ds["u"] ** 2 + ds["v"] ** 2) ** 0.5

        if half_precision:
            speeds = speeds.astype(np.float32)

        return _splrep(speeds, "height")

    def _estimate_dataset(
        self,
        height: int,
        years: Optional[slice] = None,
        months: Optional[slice] = None,
        xs: Optional[slice] = None,
        ys: Optional[slice] = None,
        use_real_data: Optional[bool] = False,
    ) -> xr.Dataset:
        params = xr.open_mfdataset(self.files).transpose("height", ...)

        if not (xs is None or ys is None):
            params = params.sel(latitude=ys, longitude=xs)

        if not (years is None or months is None):
            if months is None:
                months = slice(1, 13)
            params = params.sel(
                valid_time=get_daterange(years, months),
            )

        # If the height is in the list of known heights, we can directly return
        # the wind speed to save computation time.
        if float(height) in LEVEL_TO_HEIGHT.values():
            params = params.sel(height=height)
            return (
                ((params["u"] ** 2 + params["v"] ** 2) ** 0.5)
                .drop("height")
                .drop("model_level")
            )

        params = params[["c"]]
        spline_params = sinterp.BSpline(
            params.attrs.get("t"), params.get("c").values, k=3
        )

        params = params.drop_dims("height")
        return xr.DataArray(
            spline_params(height), dims=params.dims, coords=params.coords
        )

    def _estimate_cutout(
        self,
        height: int,
        years: slice,
        months: Optional[slice] = None,
        xs: Optional[slice] = None,
        ys: Optional[slice] = None,
        use_real_data: Optional[bool] = False,
    ) -> xr.Dataset:
        params = xr.open_mfdataset(self.files).transpose("height", ...)

        if not (xs is None or ys is None):
            params = params.sel(latitude=ys, longitude=xs)

        if not (years is None or months is None):
            if months is None:
                months = slice(1, 13)
            params = params.sel(
                valid_time=get_daterange(years, months),
            )

        if float(height) in LEVEL_TO_HEIGHT.values() and use_real_data:
            params = params.sel(height=height)
            return (
                ((params["u"] ** 2 + params["v"] ** 2) ** 0.5)
                .drop("height")
                .drop("model_level")
            )

        params = params[["c"]]
        spline_params = sinterp.BSpline(
            params.attrs.get("t"), params.get("c").values, k=3
        )

        params = params.drop_dims("height")
        return xr.DataArray(
            spline_params(height), dims=params.dims, coords=params.coords
        )
