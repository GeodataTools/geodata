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

import abc
import dataclasses
import hashlib
import itertools
import logging
from collections.abc import Sequence
from pathlib import Path

import pandas as pd
import xarray as xr
from tqdm.auto import tqdm

from ..config import DATASET_ROOT_PATH
from ..types import BoundRange, DateRange

logger = logging.getLogger(__name__)


@dataclasses.dataclass
class AtomicDataset:
    """AtomicDataset is a class that encapsulates an individual xarray file that was
    downloaded. It provides a streamlined workflow for downloading, preprocessing,
    and integrity checking of these datasets.
    """

    year: int
    month: int
    dataset: "BaseDataset"
    day: int | None = None
    file_hash: str | None = None
    url: str | None = None
    spinup: bool | None = None

    def __post_init__(self):
        if not isinstance(self.dataset, BaseDataset):
            raise ValueError("dataset must be an instance of BaseDataset")

    @property
    def path(self):
        """The path where the file should be saved"""
        if self.day is None:
            return self.dataset.storage_root / str(self.year) / f"{self.month:02d}.nc"
        else:
            return (
                self.dataset.storage_root
                / str(self.year)
                / f"{self.month:02d}"
                / f"{self.day:02d}.nc"
            )

    def check(self, integrity: bool = True):
        """Check the presence of the file and its integrity.

        Args:
            integrity: A boolean flag indicating whether to check the integrity
            of the file. If True, the file will be checked against its hash.
            If False, only the presence of the file will be checked.

        Returns:
            True if the file is present and its integrity is intact, False otherwise.
        """

        if not self.path.exists():
            logger.debug(f"{self.path} does not exist")
            return False

        if not integrity:
            return True

        # In case the file just got downloaded
        if self.file_hash is None:
            self.file_hash = self._compute_hash()
            return True

        if self.file_hash != self._compute_hash():
            logger.warning(f"{self.path} is corrupted")
            return False

        return True

    def _compute_hash(self):
        """Compute the hash of the file.

        Returns:
            The hash of the file.
        """

        hash_func = hashlib.sha256()
        with open(self.path, "rb") as f:
            while chunk := f.read(8192):  # Read file in chunks
                hash_func.update(chunk)
        return hash_func.hexdigest()


class BaseDataset(abc.ABC):
    """Dataset is a class that encapsulates any datasets natively supported
    by geodata. It provides a streamlined workflow for downloading, preprocessing,
    and storing of these datasets.

    Args:
        years: A slice object or sequence of two integers representing the
            range of years to download.
        months: A slice object or sequence of two integers representing the
            range of months to download.
        bounds (optional): A tuple of four floats representing the bounding box
            (lon_min, lat_min, lon_max, lat_max)
        **kwargs: Additional keyword arguments that are passed to the dataset.

    Notes:
        - Subclasses of BaseDataset must define the following attributes:
            - module: The module of the dataset.
            - weather_config: The configuration of the dataset.
        - Subclasses of BaseDataset must also implement the following methods:
            - download: Method to download the dataset.
            - _extra_setup: Method to handle any extra setup that is required
                for the dataset.
        - By default, the files downloaded by the dataset are defined to by monthly in
            nature. That is,  Subclasses can override this behavior by setting the `frequency`

    """

    module: str
    weather_config: str
    frequency: str = "monthly"

    def __init__(
        self,
        years: DateRange,
        months: DateRange,
        bounds: BoundRange | None = None,
        **kwargs,
    ):
        if not hasattr(self, "module"):
            raise ValueError("Subclasses of BaseDataset must define a module attribute")

        if not hasattr(self, "weather_config"):
            raise ValueError(
                "Subclasses of BaseDataset must define a weather_config attribute"
            )

        if not isinstance(years, slice):
            if isinstance(years, Sequence):
                if not all(isinstance(year, int) for year in years):
                    raise ValueError("years must be a sequence of integers")
                elif not len(years) == 2:
                    raise ValueError("years must be a sequence of length 2")
                years = slice(years[0], years[1])
            else:
                raise ValueError(
                    f"""Invalid input {years} for years. Years must either be a
                    sequence of integers or a slice object"""
                )
        self.years = years

        if not isinstance(months, slice):
            if isinstance(months, Sequence):
                if not all(isinstance(month, int) for month in months):
                    raise ValueError("months must be a sequence of integers")
                elif not len(months) == 2:
                    raise ValueError("months must be a sequence of length 2")
                elif not all(1 <= month <= 12 for month in months):
                    raise ValueError(
                        "months must be a sequence of integers between 1 and 12"
                    )
                months = slice(months[0], months[1])
            else:
                raise ValueError(
                    f"""Invalid input {months} for months. Months must either be a
                    sequence of integers or a slice object"""
                )
        self.months = months

        if bounds is not None:
            if not all(isinstance(bound, (int, float)) for bound in bounds):
                raise ValueError("bounds must be a sequence of integers or floats")
            if not len(bounds) == 4:
                raise ValueError("bounds must be a sequence of length 4")
            if not all(-180 <= bound <= 180 for bound in [bounds[0], bounds[2]]):
                raise ValueError("Longitude bounds must be between -180 and 180")
            if not all(-90 <= bound <= 90 for bound in [bounds[1], bounds[3]]):
                raise ValueError("Latitude bounds must be between -90 and 90")
        self.bounds = bounds

        self.storage_root = (
            Path(kwargs.get("dataset_root", DATASET_ROOT_PATH))
            / self.module
            / self.weather_config
        )
        if not self.storage_root.exists():
            logger.info(
                f"Storage directory for {self.__class__.__name__} does not exist, "
                f"creating now at {self.storage_root}"
            )
            self.storage_root.mkdir(parents=True)

        self._extra_kwargs = kwargs
        self._extra_setup(**kwargs)

    def _extra_setup(self, **kwargs):
        """Method to be implemented by subclasses to handle any extra setup
        that is required for the dataset.
        """

    @staticmethod
    def apply(func: callable, *args, **kwargs):
        """Method to apply a function to each file in the dataset.

        Args:
            func: A function that takes an xarray.Dataset as its first argument.
            *args: Additional arguments to pass to the function.
            **kwargs: Additional keyword arguments to pass to the function.
        """

        def wrapper(self, *args, **kwargs):
            for file in self.catalog:
                with xr.open_dataset(file.path, chunks="auto") as ds:
                    ds = func(ds, *args, **kwargs)
                    ds.to_netcdf(file.path)

        return wrapper

    def _generate_manifest(self):
        """Generate a manifest file for the dataset. This file contains
        metadata about the dataset, including the file paths and their
        integrity checks.

        TODO! This method is not ready yet!

        Returns:
            A list of dictionaries containing the metadata of each file in the
            dataset.
        """

        manifest = []
        for file in self.catalog:
            manifest.append(
                {
                    "path": str(file["save_path"]),
                    "integrity": file["downloaded"],
                }
            )

        return manifest

    @property
    def downloaded(self):
        """A boolean flag indicating whether the dataset has been prepared
        for use. This typically means that the dataset has been downloaded,
        preprocessed, and stored in a format that is ready for use.

        The basic implementation of this method checks the presence of each file in
        the catalog in the storage root. Subclasses can override this behavior if
        a more comprehensive check is required.
        """

        return all((file.check() for file in self.catalog))

    @abc.abstractmethod
    def _download_file(self, file: dict):
        """Method to download a single file from the dataset. This method
        should download the file and save it to the appropriate location.

        Args:
            file: A dictionary containing the metadata of the file to download. At the
            minimum, this dictionary should contain the following keys:
                - year: the year of the file
                - month: the month of the file
                - day: the day of the file (if applicable)
                - hour: the hour of the file (if applicable)
                - save_path: the path where the file should be saved
        """

    def download(self, force: bool = False):
        """Method to download the dataset. This method should download the
        dataset files and store them in the appropriate location.

        Args:
            force: A boolean flag indicating whether to force the download of
                the dataset, even if it has already been downloaded.
        """

        if self.downloaded and not force:
            logger.info(f"{self} has already been downloaded.")
            return

        for file in tqdm(self.catalog, unit="file", dynamic_ncols=True):
            self._download_file(file)

        logger.info(f"Downloaded {self}")
        logger.info("Cleaning and renaming coordinates")

        self.apply(self._dataset_postprocess)
        self.apply(self._rename_and_clean_coords)

    def _dataset_postprocess(self, ds: xr.Dataset | xr.DataArray, **kwargs):
        """Method to postprocess the dataset after it has been downloaded.
        This method should be implemented by subclasses to handle any
        additional processing that is required for the dataset.

        Args:
            ds: The dataset to postprocess.
            **kwargs: Additional keyword arguments to pass to the function.
        """

        return ds

    @apply
    def trim_variables(
        self, variables: Sequence[str] | None = None, **kwargs
    ) -> xr.Dataset | xr.DataArray:
        """Method to trim the dataset to only include the specified variables.

        Args:
            variables: A sequence of strings representing the variables to keep.
                If None, we will keep the variables specified in the `variables`
                attribute of the dataset.
        """

        if variables is None:
            if not hasattr(self, "variables"):
                raise ValueError(
                    "The dataset does not have a `variables` attribute defined."
                    "Please specify the variables to keep."
                )
            variables: Sequence[str] = getattr(self, "variables")

        return kwargs["ds"][variables]

    def __repr__(self):
        return "<Dataset Module={} Config={} Years={}-{} Months={}-{} {}{}>".format(
            self.module,
            self.weather_config,
            self.years.start,
            self.years.stop,
            self.months.start,
            self.months.stop,
            "Downloaded" if self.downloaded else "Not Downloaded",
            " " + self.extra_repr if self.extra_repr else "",
        )

    @property
    @abc.abstractmethod
    def projection(self):
        """The projection of the dataset. This should be a string that
        represents the projection of the dataset.
        """

    @property
    def extra_repr(self):
        return ""

    @property
    def testing(self):
        """A boolean flag indicating whether the dataset is being used for
        testing. Under this mode, only the first few days or months of the dataset
        will be downloaded (depending on the granularity). This is useful for
        testing the dataset without downloading the entire dataset.
        """

        if "testing" not in self._extra_kwargs:
            return False
        if not isinstance(self._extra_kwargs["testing"], bool):
            raise ValueError("testing must be a boolean flag")
        return self._extra_kwargs["testing"]

    @property
    def submodule(self):
        """The submodule of the dataset. This can be defined by the dataset
        using the `weather_config` attribute. If not defined, it will default
        to the name of the dataset class.
        """
        return getattr(self, "weather_config", self.__class__.__name__)

    @property
    def catalog(self) -> list["AtomicDataset"]:
        """A generator that yields all the files that need to be downloaded.
        Each iteration should return a dictionary with the following keys:
            - year: the year of the file
            - month: the month of the file
            - day: the day of the file (if applicable)
            - hour: the hour of the file (if applicable)
            - save_path: the path where the file should be saved
        """

        match self.frequency:
            case "monthly":
                cat = self._monthly_catalog()
            case "daily":
                cat = self._daily_catalog()
            case "hourly":
                cat = self._hourly_catalog()
            case _:
                raise ValueError(
                    f"Invalid frequency {self.frequency} defined for this dataset."
                )

        return cat

    def _monthly_catalog(self):
        catalog = []

        for year, month in itertools.product(
            range(self.years.start, self.years.stop + 1),
            range(self.months.start, self.months.stop + 1),
        ):
            catalog.append(AtomicDataset(year, month, dataset=self))

        return catalog

    def _daily_catalog(self) -> list["AtomicDataset"]:
        catalog = []

        for year, month in itertools.product(
            range(self.years.start, self.years.stop + 1),
            range(self.months.start, self.months.stop + 1),
        ):
            for day in range(1, pd.Timestamp(f"{year}-{month}-1").days_in_month + 1):
                catalog.append(AtomicDataset(year, month, day, dataset=self))

        return catalog

    # def _hourly_catalog(self):
    #     catalog = []

    #     for year, month in itertools.product(
    #         range(self.years.start, self.years.stop + 1),
    #         range(self.months.start, self.months.stop + 1),
    #     ):
    #         for day in range(1, pd.Timestamp(f"{year}-{month}-1").days_in_month + 1):
    #             for hour in range(24):
    #                 save_path = (
    #                     self.storage_root
    #                     / f"{year}_{month:02d}_{day:02d}_{hour:02d}.nc"
    #                 )
    #                 catalog.append(
    #                     {
    #                         "year": year,
    #                         "month": month,
    #                         "day": day,
    #                         "hour": hour,
    #                         "save_path": save_path,
    #                     }
    #                 )
    #     return catalog

    def _rename_and_clean_coords(ds: xr.Dataset, add_lon_lat: bool = True):
        """Rename 'lon'/'longitude' and 'lat'/'latitude' columns to 'x' and 'y'

        Optionally (add_lon_lat, default:True) preserves latitude and longitude columns as 'lat' and 'lon'.

        Args:
            ds (xarray.Dataset): Dataset to rename
            add_lon_lat (bool, optional): Add lon/lat columns. Defaults to True.

        Returns:
            xarray.Dataset: Dataset with renamed coordinates
        """

        # Rename latitude / lat -> y, longitude / lon -> x
        if "latitude" in list(ds.coords):
            ds = ds.rename({"latitude": "y"})
        if "longitude" in list(ds.coords):
            ds = ds.rename({"longitude": "x"})
        if "lat" in list(ds.coords):
            ds = ds.rename({"lat": "y"})
        if "lon" in list(ds.coords):
            ds = ds.rename({"lon": "x"})

        if add_lon_lat:
            ds = ds.assign_coords(lon=ds.coords["x"], lat=ds.coords["y"])

        return ds
