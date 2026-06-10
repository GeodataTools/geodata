# Copyright 2024-2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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

import logging
from typing import Hashable
from typing import cast

import numpy as np
import scipy.interpolate as sinterp
import xarray as xr
from xarray.namedarray.pycompat import array_type

from ...logging import logger
from ...utils import rechunk_dataset
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


def _make_interp_coeff(x, y, k, t, check_finite=False):
    """Compute interpolation coefficients for a single block.
    
    Args:
        x: 1D array of x coordinates (must be C-contiguous)
        y: Data array to interpolate
        k: Spline degree
        t: Knot vector
        check_finite: Whether to check for finite values
        
    Returns:
        Spline coefficients
    """
    logger.debug(f"[_make_interp_coeff] Called with x type: {type(x)}, x shape: {np.asarray(x).shape if hasattr(x, 'shape') else 'no shape'}, "
                 f"x dtype: {np.asarray(x).dtype if hasattr(x, 'dtype') else type(x)}, "
                 f"y type: {type(y)}, y shape: {np.asarray(y).shape if hasattr(y, 'shape') else 'no shape'}, "
                 f"k: {k}, t type: {type(t)}, t shape: {np.asarray(t).shape if hasattr(t, 'shape') else 'no shape'}")
    
    # Ensure x is C-contiguous and writable
    x = _memoryview_safe(np.asarray(x, dtype=float))
    logger.debug(f"[_make_interp_coeff] After _memoryview_safe: x shape: {x.shape}, x dtype: {x.dtype}, x.flags.c_contiguous: {x.flags.c_contiguous}")
    
    try:
        result = sinterp.make_interp_spline(x, y, k=k, t=t, check_finite=check_finite).c
        logger.debug(f"[_make_interp_coeff] Successfully computed coefficients, shape: {result.shape}")
        return result
    except Exception as e:
        logger.error(f"[_make_interp_coeff] ERROR in make_interp_spline: {type(e).__name__}: {e}")
        logger.error(f"[_make_interp_coeff] x details: shape={x.shape}, dtype={x.dtype}, c_contiguous={x.flags.c_contiguous}")
        logger.error(f"[_make_interp_coeff] y details: type={type(y)}, shape={np.asarray(y).shape if hasattr(y, 'shape') else 'N/A'}")
        logger.error(f"[_make_interp_coeff] t details: type={type(t)}, shape={np.asarray(t).shape if hasattr(t, 'shape') else 'N/A'}")
        raise


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

    logger.debug(f"[_splrep] Starting with dim={dim}, k={k}, a shape: {a.shape}, a dims: {a.dims}")
    
    # Make sure that dim is on axis 0
    a = a.transpose(dim, ...)
    x: np.ndarray = a.coords[dim].values
    logger.debug(f"[_splrep] After transpose: a shape: {a.shape}, x shape: {x.shape}, x dtype: {x.dtype}")

    if x.dtype.kind == "M":
        # Same treatment will be applied to x_new.
        # Allow x_new.dtype==M8[D] and x.dtype==M8[ns], or vice versa
        x = x.astype("M8[ns]").astype(float)
        logger.debug(f"[_splrep] Converted datetime x to float, new dtype: {x.dtype}")
    
    # Ensure x is C-contiguous and properly typed
    x = _memoryview_safe(np.asarray(x, dtype=float))
    logger.debug(f"[_splrep] After _memoryview_safe: x shape: {x.shape}, x dtype: {x.dtype}, x.flags.c_contiguous: {x.flags.c_contiguous}")

    t = sinterp._bsplines._not_a_knot(x, k=k)
    logger.debug(f"[_splrep] Computed knots t, shape: {t.shape}, dtype: {t.dtype}")

    if isinstance(a.data, array_type("dask")):
        from dask.array import map_blocks
        from dask.diagnostics.progress import ProgressBar

        logger.debug(f"[_splrep] Data is dask array, chunks: {a.data.chunks}, shape: {a.data.shape}")

        if len(a.data.chunks[0]) > 1:
            logger.debug(f"[_splrep] Rechunking dimension {dim} to -1 (was: {a.data.chunks[0]})")
            a = a.chunk({dim: -1})
            logger.debug(f"[_splrep] After rechunking, chunks: {a.data.chunks}")

        pbar = ProgressBar()
        if logger.level <= logging.INFO:
            pbar.register()

        # Create a wrapper function that captures x and t as closures
        # This ensures they're passed correctly to each block
        def _block_interp_coeff(y_block, x=x, k=k, t=t, check_finite=False):
            y_block = np.asarray(y_block)
            logger.debug(f"[_block_interp_coeff] Called with y_block type: {type(y_block)}, y_block shape: {y_block.shape}, "
                        f"x type: {type(x)}, x shape: {x.shape if hasattr(x, 'shape') else 'no shape'}, "
                        f"t type: {type(t)}, t shape: {t.shape if hasattr(t, 'shape') else 'no shape'}")
            
            # Handle empty blocks - return empty array with correct shape
            if y_block.size == 0 or any(s == 0 for s in y_block.shape):
                logger.debug(f"[_block_interp_coeff] Empty block detected, returning empty array with shape: {y_block.shape}")
                # Return empty array with same shape as input (coefficients have same shape as input)
                return np.empty_like(y_block, dtype=float)
            
            try:
                result = _make_interp_coeff(x, y_block, k=k, t=t, check_finite=check_finite)
                logger.debug(f"[_block_interp_coeff] Successfully computed, result shape: {result.shape}")
                return result
            except Exception as e:
                logger.error(f"[_block_interp_coeff] ERROR: {type(e).__name__}: {e}")
                raise

        logger.debug(f"[_splrep] Calling map_blocks with a.data shape: {a.data.shape}, chunks: {a.data.chunks}")
        logger.debug(f"[_splrep] x closure value: shape={x.shape}, dtype={x.dtype}, c_contiguous={x.flags.c_contiguous}")
        logger.debug(f"[_splrep] t closure value: shape={t.shape}, dtype={t.dtype}")
        
        try:
            c = map_blocks(
                _block_interp_coeff,
                a.data,
                dtype=float,
                drop_axis=[],
            )
            logger.debug(f"[_splrep] map_blocks returned, c type: {type(c)}, c shape: {c.shape if hasattr(c, 'shape') else 'N/A'}")
        except Exception as e:
            logger.error(f"[_splrep] ERROR in map_blocks: {type(e).__name__}: {e}")
            raise

        if logger.level <= logging.INFO:
            pbar.unregister()

    else:
        logger.debug(f"[_splrep] Data is numpy array (not dask), shape: {a.data.shape}, dtype: {a.data.dtype}")
        c = _make_interp_coeff(x, a.data, k=k, t=t, check_finite=False)
        logger.debug(f"[_splrep] Computed coefficients (numpy), shape: {c.shape}")

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


def _splev_ker(c: np.ndarray, t: np.ndarray, k: int, height: np.ndarray) -> np.ndarray:
    return np.atleast_1d(sinterp.splev(height, (t, c, k)))


def _splev(da: xr.Dataset, height: float) -> xr.DataArray:
    height_arr = np.atleast_1d(height)
    return xr.apply_ufunc(
        _splev_ker,
        da["c"],
        input_core_dims=[["height"]],
        output_core_dims=[[]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[da["c"].dtype],
        kwargs={"t": da.attrs["t"], "k": da.attrs["k"], "height": height_arr},
    )


class WindInterpolationModel(WindBaseModel):
    """Wind speed estimation based on a spline interpolation of the wind speed at different heights.

    This model uses the ERA5 3D dataset to estimate wind speed at a given height using
    spline interpolation.

    For ``estimate(..., xs=..., ys=...)``, each spatial slice may use either bound order
    (``slice(low, high)`` or ``slice(high, low)``); ``BaseModel.estimate`` normalizes
    slices to the coordinate monotonic direction before ``xarray.Dataset.sel``, including
    for descending coordinates (e.g. latitude).

    Example:

    >>> from geodata import Dataset
    >>> from geodata.model.wind import WindInterpolationModel
    >>> dataset = Dataset(module="era5", weather_data_config="wind_3d_hourly", years=slice(2010, 2010), months=slice(1, 2))
    >>> model = WindExtrapolationModel(dataset)
    >>> model.prepare()
    >>> model.estimate(height=12, xs=slice(1, 2), ys=slice(1, 2), years=slice(2010, 2010), months=slice(1, 2))
    """

    SUPPORTED_WEATHER_DATA_CONFIGS = ("wind_3d_hourly", "wind_3d_hourly_test")

    def _prepare_dataset(
        self,
        source: xr.Dataset,
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
            "model_level" in source.coords
        ), "Dataset does not contain model levels. Please double-check the dataset."

        source.coords["model_level"] = np.array(
            [LEVEL_TO_HEIGHT[int(level)] for level in source["model_level"].values]
        )
        source = (
            source.rename({"model_level": "height"})
            .transpose("height", ...)
            .sortby("height")
        )

        logger.debug("Shape of heights: %s", source["height"].shape)
        speeds = (source["u"] ** 2 + source["v"] ** 2) ** 0.5
        logger.debug(f"[_prepare_dataset] Computed speeds, shape: {speeds.shape}, dims: {speeds.dims}, "
                    f"is dask: {isinstance(speeds.data, array_type('dask'))}")

        logger.info(f"[_prepare_dataset] Calling _splrep with speeds shape: {speeds.shape}")
        params = _splrep(speeds, "height")
        logger.info(f"[_prepare_dataset] _splrep returned params, type: {type(params)}, data_vars: {list(params.data_vars.keys())}")
        if half_precision:
            params = params.astype(np.float32)

        return params

    def _estimate_dataset(self, params: xr.Dataset, **kwargs) -> xr.DataArray:
        height = float(kwargs["height"])
        params = params.transpose("height", ...)
        params = cast(
            xr.Dataset,
            rechunk_dataset(params, force_full_chunk_dims=["height"]),
        )
        result = _splev(params, height)

        # Some upstream ERA5 pipelines historically use `valid_time` for the time-like
        # coordinate. Normalize to `time` so we can enforce consistent dims.
        if "valid_time" in result.dims or "valid_time" in result.coords:
            result = result.rename({"valid_time": "time"})

        # Standardize output dimension order across wind models:
        # `("time", "x", "y")` (wind interpolation should match pvlib).
        desired_order = ("time", "x", "y")
        if all(d in result.dims for d in desired_order):
            result = result.transpose(*desired_order)
        return result
