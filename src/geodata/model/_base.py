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
import shutil
from typing import Optional

import xarray as xr
from tqdm.auto import tqdm

from ..config import model_dir
from ..datasets._base import BaseDataset
from ..logging import logger
from .results import DailyModelResult, MonthlyModelResult, ResultType

if importlib.util.find_spec("h5netcdf") is not None:
    XR_PARALLEL = True
    XR_ENGINE = "h5netcdf"
else:
    XR_PARALLEL = False
    XR_ENGINE = None
    logger.warning(
        "h5netcdf is not installed. Parallel reading of netCDF files will be disabled. "
        "This could have some performance implications."
    )

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

    SUPPORTED_WEATHER_DATA_CONFIGS: tuple[str]

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
        params = xr.open_mfdataset(files, engine=XR_ENGINE, parallel=XR_PARALLEL)

        if xs is not None:
            params = params.sel(x=xs)
        if ys is not None:
            params = params.sel(y=ys)

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

                with xr.open_mfdataset(
                    result.ref_files, engine=XR_ENGINE, parallel=XR_PARALLEL
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
