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

import abc
import itertools
import logging
from collections.abc import Sequence
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm

from ..config import DATASET_ROOT_PATH

logger = logging.getLogger(__name__)


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
        years: Sequence[int] | slice,
        months: Sequence[int] | slice,
        bounds: Sequence[int] | None = None,
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
                f"""Storage directory for {self.__class__.__name__}
                does not exist, creating now at {self.storage_root}"""
            )
            self.storage_root.mkdir(parents=True)

        self._extra_setup(**kwargs)

    @abc.abstractmethod
    def _extra_setup(self, **kwargs):
        """Method to be implemented by subclasses to handle any extra setup
        that is required for the dataset.
        """

    @property
    def downloaded(self):
        """A boolean flag indicating whether the dataset has been prepared
        for use. This typically means that the dataset has been downloaded,
        preprocessed, and stored in a format that is ready for use.

        The basic implementation of this method checks the presence of each file in
        the catalog in the storage root. Subclasses can override this behavior if
        a more comprehensive check is required.
        """

        return all((file["save_path"].exists() for file in self.catalog))

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

        for file in tqdm(
            self.catalog, desc="Downloading", unit="file", dynamic_ncols=True
        ):
            self._download_file(file)

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
    def extra_repr(self):
        return ""

    @property
    def submodule(self):
        """The submodule of the dataset. This can be defined by the dataset
        using the `weather_config` attribute. If not defined, it will default
        to the name of the dataset class.
        """
        return getattr(self, "weather_config", self.__class__.__name__)

    @property
    def catalog(self):
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
                return self._monthly_catalog()
            case "daily":
                return self._daily_catalog()
            case "hourly":
                return self._hourly_catalog()
            case _:
                raise ValueError(
                    f"Invalid frequency {self.frequency} defined for this dataset."
                )

    def _monthly_catalog(self):
        catalog = []

        for year, month in itertools.product(
            range(self.years.start, self.years.stop + 1),
            range(self.months.start, self.months.stop + 1),
        ):
            save_path = self.storage_root / f"{year}_{month:02d}.nc"
            catalog.append({"year": year, "month": month, "save_path": save_path})

        return catalog

    def _daily_catalog(self):
        catalog = []

        for year, month in itertools.product(
            range(self.years.start, self.years.stop + 1),
            range(self.months.start, self.months.stop + 1),
        ):
            for day in range(1, pd.Timestamp(f"{year}-{month}-1").days_in_month + 1):
                save_path = self.storage_root / f"{year}_{month:02d}_{day:02d}.nc"
                catalog.append(
                    {"year": year, "month": month, "day": day, "save_path": save_path}
                )

        return catalog

    def _hourly_catalog(self):
        catalog = []

        for year, month in itertools.product(
            range(self.years.start, self.years.stop + 1),
            range(self.months.start, self.months.stop + 1),
        ):
            for day in range(1, pd.Timestamp(f"{year}-{month}-1").days_in_month + 1):
                for hour in range(24):
                    save_path = (
                        self.storage_root
                        / f"{year}_{month:02d}_{day:02d}_{hour:02d}.nc"
                    )
                    catalog.append(
                        {
                            "year": year,
                            "month": month,
                            "day": day,
                            "hour": hour,
                            "save_path": save_path,
                        }
                    )
        return catalog
