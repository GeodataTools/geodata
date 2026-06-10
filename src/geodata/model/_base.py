# Copyright 2023-2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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


import abc
import importlib.util
import os
import platform
import shutil
from collections.abc import Collection
from typing import ClassVar, Optional

import numpy as np
import xarray as xr
from tqdm.auto import tqdm

from ..config import model_dir
from ..datasets._base import BaseDataset
from ..logging import logger
from .results import DailyModelResult, MonthlyModelResult, ResultType

if importlib.util.find_spec("h5netcdf") is not None:
    XR_ENGINE = "h5netcdf"
    XR_PARALLEL_DEFAULT = True
else:
    XR_PARALLEL_DEFAULT = False
    XR_ENGINE = None
    logger.warning(
        "h5netcdf is not installed. Parallel reading of netCDF files will be disabled. "
        "This could have some performance implications."
    )


def _normalize_slice_for_sel(coord: xr.DataArray, s: slice) -> slice:
    """Return a slice for ``.sel()`` that matches the coordinate direction.

    xarray's ``.sel(dim=slice(a, b))`` returns empty when the dimension is descending
    (e.g. ERA5 latitude) or when the user passes ``slice(high, low)`` on an ascending
    dimension. This helper interprets the slice as the inclusive logical range
    ``[min(start, stop), max(start, stop)]`` and returns bounds in the order required
    by ``.sel()`` for that coordinate's monotonic direction.
    """
    if not isinstance(s, slice) or s.step not in (None, 1):
        return s
    if s.start is None or s.stop is None:
        return s
    lo, hi = min(s.start, s.stop), max(s.start, s.stop)
    vals = np.asarray(coord.values).ravel()
    if len(vals) < 2:
        return slice(lo, hi)
    descending = np.all(np.diff(vals) <= 0)
    if descending:
        return slice(hi, lo)
    return slice(lo, hi)


def _is_in_dask_worker_on_linux() -> bool:
    """Check if we're running in a Dask worker process on Linux.
    
    Returns:
        bool: True if we're in a Dask worker on Linux, False otherwise.
    """
    if platform.system() != "Linux":
        return False
    
    try:
        from dask.distributed import get_worker
        try:
            get_worker()
            return True
        except ValueError:
            # Not in a worker process
            return False
    except ImportError:
        # dask.distributed not available
        return False


def _is_dask_using_processes_on_linux() -> bool:
    """Check if Dask is being used with processes on Linux.
    
    Returns:
        bool: True if Dask is using processes on Linux, False otherwise.
        
    Note:
        This checks if there's an active Dask client using processes.
        When Dask uses processes, h5netcdf has issues with HDF5 dimension scales.
    """
    if platform.system() != "Linux":
        return False
    
    try:
        from dask.distributed import get_client, get_worker
        try:
            get_client()
            # Check if we're in a worker (which means processes are being used)
            try:
                get_worker()
                return True
            except ValueError:
                # Not in a worker, but check if client exists and might use processes
                # We can't easily detect this from the main process, so we'll be conservative
                # and assume processes might be used if a client exists
                # The actual check will happen in workers via _is_in_dask_worker_on_linux
                return False
        except ValueError:
            # No active client
            return False
    except ImportError:
        # dask.distributed not available
        return False


def _get_xr_engine() -> str | None:
    """Get the appropriate xarray engine to use for opening NetCDF files.
    
    Returns:
        str | None: The engine name to use, or None for default.
    """
    logger.debug(f"_get_xr_engine: Returning engine {XR_ENGINE}")
    return XR_ENGINE


def _should_use_parallel_reading() -> bool:
    """Determine if parallel reading should be used for xarray open_mfdataset.
    
    Returns:
        bool: True if parallel reading should be used, False otherwise.
        
    Note:
        Parallel reading is disabled when Dask is using processes on Linux,
        as h5netcdf has issues with HDF5 dimension scales in that case.
    """
    if not XR_PARALLEL_DEFAULT:
        logger.debug("_should_use_parallel_reading: XR_PARALLEL_DEFAULT is False, returning False")
        return False
    
    in_worker = _is_in_dask_worker_on_linux()
    using_processes = _is_dask_using_processes_on_linux()
    
    if in_worker or using_processes:
        logger.info(
            f"_should_use_parallel_reading: Disabling parallel reading "
            f"(in_worker={in_worker}, using_processes={using_processes})"
        )
        return False
    
    logger.debug(f"_should_use_parallel_reading: Returning {XR_PARALLEL_DEFAULT}")
    return XR_PARALLEL_DEFAULT

# Parse the MAX_WORKERS environment variable if present
MAX_WORKERS = os.getenv("MAX_WORKERS")
if MAX_WORKERS is not None:
    try:
        max_workers = int(MAX_WORKERS)
    except ValueError:
        logger.warning(
            "MAX_WORKERS environment variable is not an integer. Using default value."
        )
        MAX_WORKERS = None


class BaseModel(abc.ABC):
    """Base class for geospatial modeling.

    Args:
        name (str): The name of the model.
        source (BaseDataset): The source of the model.
        interpolate (bool, optional): Interpolate the source to the same grid as the target. Defaults to False.
        quick_check (bool, optional): Quick check for the model. Defaults to False. If True, the model parameters will be checked for presence, but not the integrity.
        **kwargs: Additional keyword arguments to pass to the model.
    """

    SUPPORTED_WEATHER_DATA_CONFIGS: ClassVar[Collection[str]]

    def __init__(self, source: BaseDataset, **kwargs):
        if source.weather_config not in self.SUPPORTED_WEATHER_DATA_CONFIGS:
            raise ValueError(
                f"Weather data config {source.weather_config} is not supported by this model."
            )

        if not source.downloaded:
            raise ValueError("The source Dataset for this model is not prepared.")

        self.source = source
        self.quick_check = kwargs.get("quick_check", False)
        self._extra_kwargs = kwargs
        self._prepared = False

        self._ref_path = model_dir.parent / self.source.module
        self._results: dict[int, dict[int, ResultType]] = self._prepare_results()

    def __repr__(self):
        return f"Model(source={self.source}, type={self.type})"

    @property
    def frequency(self) -> str:
        """Frequency of the model."""
        return self.source.frequency

    @property
    @abc.abstractmethod
    def type(self) -> str:
        """Type of the model."""

    def _prepare_results(self) -> dict[int, dict[int, ResultType]]:
        """Prepare the results of the model.

        Returns:
            dict: Dictionary with the results of the model.
        """

        years = list(range(self.source.years.start, self.source.years.stop + 1))
        months = list(range(self.source.months.start, self.source.months.stop + 1))

        results: dict[int, dict[int, ResultType]] = {}
        for year in years:
            results[year] = {}
            for month in months:
                match self.frequency:
                    case "daily" | "hourly":
                        results[year][month] = DailyModelResult.from_year_month(
                            self, year, month
                        )
                    case "monthly":
                        results[year][month] = MonthlyModelResult.from_year_month(
                            self, year, month
                        )
                results[year][month].path.mkdir(parents=True, exist_ok=True)

        return results

    @property
    def results(self):
        """Get the results of the model.

        Returns:
            dict: Dictionary with the results of the model.
        """
        return self._results

    def get_result_year_month(self, years: slice, months: slice) -> list[ResultType]:
        """Get the result of the model for a given year and month range.

        Args:
            years (slice): Year range.
            months (slice): Month range.
        Returns:
            list: List of DailyModelResult objects.
        """
        year_d = [self._results[y] for y in range(years.start, years.stop + 1)]
        return [y[m] for y in year_d for m in range(months.start, months.stop + 1)]

    @property
    def flattened_results(self) -> list[ResultType]:
        """Flatten the results of the model.

        Returns:
            list: List of ModelResult objects.
        """

        return [
            self._results[year][month]
            for year in self._results
            for month in self._results[year]
        ]

    def estimate(
        self,
        years: Optional[slice] = None,
        months: Optional[slice] = None,
        xs: Optional[slice] = None,
        ys: Optional[slice] = None,
        **kwargs,
    ) -> xr.DataArray:
        """Estimate the wind speed at given coordinates.

        Args:
            years (slice, optional): Years.
            months (slice, optional): Months. If None, all months are estimated.
            xs (slice, optional): X coordinates. If None, all x coordinates in source are estimated.
            ys (slice, optional): Y coordinates. If None, all y coordinates in source are estimated.
            **kwargs: Additional keyword arguments to pass to the model.

        Returns:
            xr.DataArray: Dataset with wind speed.
        """
        if not self.prepared:
            raise RuntimeError(
                "The model is not prepared. Please prepare the model first."
            )

        if years is None and months is None:
            results = self.flattened_results
        elif months is None:
            results = self.get_result_year_month(years, slice(1, 13))
        else:
            results = self.get_result_year_month(years, months)

        files = sum([result.files for result in results], [])
        engine = _get_xr_engine()
        parallel = _should_use_parallel_reading()
        logger.info(
            f"estimate: Opening {len(files)} files with engine={engine}, parallel={parallel}"
        )
        params = xr.open_mfdataset(files, engine=engine, parallel=parallel)

        if xs is not None:
            x_slice = (
                _normalize_slice_for_sel(params.coords["x"], xs)
                if "x" in params.coords
                else xs
            )
            params = params.sel(x=x_slice)
        if ys is not None:
            y_slice = (
                _normalize_slice_for_sel(params.coords["y"], ys)
                if "y" in params.coords
                else ys
            )
            params = params.sel(y=y_slice)

        output = self._estimate_dataset(params, **kwargs)
        params.close()
        return output

    @property
    def prepared(self) -> bool:
        """Check if the model is prepared.

        Returns:
            bool: True if prepared.
        """

        if not self._prepared:
            self._prepared = self._check_prepared()
        return self._prepared

    def _check_prepared(self) -> bool:
        for year in self._results:
            for month in self._results[year]:
                if not self._results[year][month].prepared:
                    return False
        return True

    def prepare(self, force: bool = False):
        """Prepare the model.

        Args:
            force (bool, optional): Force re-prepare the model. Defaults to False.
        """

        if self.prepared and not force:
            logger.info("The model is already prepared.")
            return

        for result in tqdm(self.flattened_results):
            if not result.prepared or force:
                shutil.rmtree(result.path, ignore_errors=True)
                result.path.mkdir(parents=True, exist_ok=True)

                engine = _get_xr_engine()
                parallel = _should_use_parallel_reading()
                logger.info(
                    f"prepare: Opening {len(result.ref_files)} files with engine={engine}, parallel={parallel}"
                )
                with xr.open_mfdataset(
                    result.ref_files,
                    engine=engine,
                    parallel=parallel,
                ) as ds:
                    prepared_ds = self._prepare_dataset(ds)
                    result.register(prepared_ds)

            result.dump()
        logger.info("Model prepared successfully.")

    @abc.abstractmethod
    def _prepare_dataset(self, source: xr.Dataset) -> xr.Dataset:
        """Prepare the parameters of a specific source dataset file.

        Args:
            source (xr.Dataset): Source dataset.

        Returns:
            xr.Dataset: Prepared parameter dataset.
        """

    @abc.abstractmethod
    def _estimate_dataset(self, params: xr.Dataset, **kwargs) -> xr.DataArray:
        """Estimate the wind speed from a dataset.

        Args:
            params (xr.Dataset): Parameters of the model.
            **kwargs: Additional keyword arguments to pass to the model.

        Returns:
            xr.DataArray: Result after modeling.
        """
