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


from typing import Hashable, Iterable, Optional, Tuple, Union

import numpy as np
import scipy.interpolate as sinterp
import xarray as xr
from dask import array as da
from tqdm.dask import TqdmCallback
from xarray.namedarray.pycompat import array_type

from ...logging import logger
from ...utils import get_daterange
from ._base import WindBaseModel

Boundary = Iterable[Tuple[int, float]]
BCType = Union[Tuple[Boundary, Boundary], str, None]


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


def _make_interp_coeffs(
    x: np.ndarray,
    y: np.ndarray,
    k: int = 3,
    t: np.ndarray | None = None,
    bc_type: BCType = None,
    axis: int = 0,
    check_finite: bool = True,
) -> np.ndarray:
    """
    Compute the knots of the B-spline.
    Borrowed from https://github.com/crusaderky/xarray_extras/blob/main/xarray_extras/kernels/interpolate.py

    Args:
        t (np.ndarray | None): Knots array, as calculated by :func:`make_interp_knots`.
            - For k=0, must always be None (the coefficients are not a function of the knots).
            - For k=1, set to None if t has been calculated by :func:`make_interp_knots`; pass a vector if it already existed before.
            - For k=2 and k=3, must always pass either the output of :func:`make_interp_knots` or a pre-generated vector.

    Returns:
        np.ndarray: Interpolation coefficients.
    """
    x = _memoryview_safe(x)
    y = _memoryview_safe(y)
    if t is not None:
        t = _memoryview_safe(t)

    return sinterp.make_interp_spline(
        x, y, k, t, bc_type=bc_type, axis=axis, check_finite=check_finite
    ).c


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

    if isinstance(a.data, array_type("dask")):
        from dask.array import map_blocks

        logger.info("Computing interpolation coefficients using Dask.")

        a = a.chunk({dim: -1})
        if len(a.data.chunks[0]) > 1:
            raise NotImplementedError(
                "Unsupported: multiple chunks on interpolation dim"
            )

        with TqdmCallback(desc="Compute Coefficients"):
            c = map_blocks(
                _make_interp_coeffs,
                x,
                a.data,
                k=k,
                t=t,
                check_finite=False,
                dtype=float,
            )

    else:
        c = _make_interp_coeffs(x, a.data, k=k, t=t, check_finite=False)

    return xr.Dataset(
        data_vars={
            "t": ("__t__", t),
            "c": (a.dims, c),
        },
        coords=a.coords,
        attrs={
            "spline_dim": dim,
            "k": k,
        },
    )


def _ker_splev(
    x_new: np.ndarray,
    t: np.ndarray,
    c: np.ndarray,
    k: int = 3,
    extrapolate: bool | str = True,
) -> np.ndarray:
    """Generate a BSpline object on the fly from knots and coefficients and
    evaluate it on x_new.

    See :class:`scipy.interpolate.BSpline` for all parameters.
    https://github.com/crusaderky/xarray_extras/blob/main/xarray_extras/kernels/interpolate.py
    """
    t = _memoryview_safe(t)
    c = _memoryview_safe(c)
    x_new = _memoryview_safe(x_new)
    spline = sinterp.BSpline.construct_fast(t, c, k, axis=0, extrapolate=extrapolate)
    return spline(x_new)


def _splev(
    x_new: object, tck: xr.Dataset, extrapolate: bool | str = True
) -> xr.DataArray:
    """
    Evaluate the B-spline generated with :func:`splrep`.
    Borrowed from https://github.com/crusaderky/xarray_extras/blob/main/xarray_extras/interpolate.py

    Args:
        x_new: Any :class:`~xr.DataArray` with any number of dims, not necessarily
            the original interpolation dim. Alternatively, it can be any 1-dimensional
            array-like; it will be automatically converted to a :class:`~xr.DataArray`
            on the interpolation dim.
        tck (xr.Dataset): As returned by :func:`splrep`. It can have been:
            - transposed (not recommended, as performance will drop if c is not C-contiguous)
            - sliced, reordered, or (re)chunked, on any dim except the interpolation dim
            - computed from dask to numpy backend
            - round-tripped to disk
        extrapolate:
            True: Extrapolate the first and last polynomial pieces of b-spline functions active on the base interval
            False: Return NaNs outside of the base interval
            'periodic': Periodic extrapolation is used
            'clip': Return y[0] and y[-1] outside of the base interval

    Returns:
        xr.DataArray: DataArray with all dims of the interpolated array, minus the interpolation dim, plus all dims of x_new

    See :func:`splrep` for usage example.
    """

    # Pre-process x_new into a DataArray
    if not isinstance(x_new, xr.DataArray):
        if not isinstance(x_new, array_type("dask")):
            x_new = np.array(x_new)
        if x_new.ndim == 0:
            dims = []
        elif x_new.ndim == 1:
            dims = [tck.spline_dim]
        else:
            raise ValueError(
                "N-dimensional x_new is only supported if x_new is a DataArray"
            )
        x_new = xr.DataArray(x_new, dims=dims, coords={tck.spline_dim: x_new})

    dim = tck.spline_dim
    t = tck.t
    c = tck.c
    k = tck.k

    invalid_dims = {*x_new.dims} & {*c.dims} - {dim}
    if invalid_dims:
        raise ValueError(
            "Overlapping dims between interpolated "
            "array and x_new: " + ",".join(str(d) for d in invalid_dims)
        )

    if t.shape != (c.sizes[dim] + k + 1,):
        raise ValueError("Interpolated dimension has been sliced")

    if x_new.dtype.kind == "M":  # datetime
        # Note that we're modifying the x_new values, not the x_new coords
        x_new = x_new.astype("M8[ns]").astype(float)
    elif x_new.dtype.kind == "m":  # timedelta
        x_new = x_new.astype("m8[ns]").astype(float)

    if extrapolate == "clip":
        x = tck.coords[dim].values

        if x.dtype.kind == "M":  # datetime
            x = x.astype("M8[ns]").astype(float)
        elif x.dtype.kind == "m":  # timedelta
            x = x.astype("m8[ns]").astype(float)

        x_new = np.clip(x_new, x[0].tolist(), x[-1].tolist())
        extrapolate = False

    if c.dims[0] != dim:
        c = c.transpose(dim, *[d for d in c.dims if d != dim])

    if any(isinstance(v.data, array_type("dask")) for v in (x_new, t, c)):
        if t.chunks and len(t.chunks[0]) > 1:
            raise NotImplementedError(
                "Unsupported: multiple chunks on interpolation dim"
            )
        if c.chunks and len(c.chunks[0]) > 1:
            raise NotImplementedError(
                "Unsupported: multiple chunks on interpolation dim"
            )

        # omitting t and c
        x_new_axes = "abdefghijklm"[: x_new.ndim]
        c_axes = "nopqrsuvwxyz"[: c.ndim - 1]

        y_new = da.blockwise(
            _ker_splev,
            x_new_axes + c_axes,
            x_new.data,
            x_new_axes,
            t.data,
            "t",
            c.data,
            "c" + c_axes,
            k=k,
            extrapolate=extrapolate,
            concatenate=True,
            dtype=float,
        )
    else:
        y_new = _ker_splev(x_new.values, t.values, c.values, k, extrapolate=extrapolate)

    y_new = xr.DataArray(y_new, dims=x_new.dims + c.dims[1:], coords=x_new.coords)
    y_new.coords.update({k: c for k, c in c.coords.items() if dim not in c.dims})

    return y_new


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

        if not (years is None and months is None):
            if months is None:
                months = slice(1, 12)
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

        return _splev(height, tck=params).astype("float32").chunk()

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

        if not (years is None and months is None):
            if months is None:
                months = slice(1, 12)
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

        return _splev(height, tck=params).astype("float32")
