# Copyright 2020 Michael Davidson (UCSD), William Honaker.
# Copyright 2020, 2023 Xiqiang Liu.

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


"""
GEODATA

Geospatial Data Collection and "Pre-Analysis" Tools
"""

import importlib
import logging
import os
import shutil
from calendar import monthrange
from collections.abc import Iterable
from tempfile import mkstemp
from typing import Optional

import numpy as np
import xarray as xr
from shapely.geometry import box

from . import config, datasets  # noqa: F401

logger = logging.getLogger(__name__)


class Dataset:
    """Dataset is a class that encapsulates any datasets natively supported by geodata.
    It provides a streamlined workflow for downloading, preprocessing, and storing of these datasets.

    Args:
        module (str): Name of dataset module
        weather_data_config (str): The type of weather data configuration to use. For more information, please
            refer to the documentation for the dataset you are using.
        years (slice): The years to download.
            For example, `slice(2000, 2010)` will download all data from 2000 to 2010.
        months (slice, optional): The months to download.
            For example, `slice(1, 12)` will download all data from January to December.
            Defaults to `slice(1, 12)`.
        opendap (bool, optional): Whether to use OpenDAP protocol for downloading.
            Defaults to `False`.
        bounds (Optional[Iterable], optional): The bounds to download.
            For example, `[-180, 90, 180, -90]` will download all data from the entire globe.
            If omitted, global data will be downloaded.
    """

    def __init__(
        self,
        module: str,
        weather_data_config: str,
        years: slice,
        months: Optional[slice] = None,
        opendap: bool = False,
        bounds: Optional[Iterable] = None,
    ):
        self.module = module
        self.config = weather_data_config
        self.dataset_module = importlib.import_module(f"geodata.datasets.{self.module}")
        self.weatherconfig = self.weather_data_config[self.config]
        self.datadir: str = self.dataset_module.datadir
        self.opendap = opendap

        self.years = years

        self.months = slice(1, 12)
        if months is None:
            logger.info("No months specified, defaulting to 1-12")
        else:
            self.months = months

        self.prepared = False
        self.toDownload = []
        self.downloadedFiles = []
        self.totalFiles = []
        incomplete_count = 0

        xs, ys = None, None
        if bounds is not None:
            self.bounds = bounds
            if isinstance(self.bounds, Iterable) and len(self.bounds) == 4:
                x1, y1, x2, y2 = bounds
                xs = slice(x1, x2)
                ys = slice(y2, y1)
            else:
                raise ValueError(
                    "Specified bounds parameter should be list with North, West, South, East coordinates."
                )
        else:
            logger.info("Bounds was not specified, default to global bounds.")
            self.bounds = None

        if os.path.isdir(self.datadir):
            # Directory for dataset exists
            logger.info("Directory %s found, checking for completeness.", self.datadir)
            self.prepared = True
            check_complete = True
        else:
            logger.info("Directory %s not found.", self.datadir)
            check_complete = False

        step = years.step if years.step else 1
        yrs = range(years.start, years.stop + step, step)
        step = months.step if months.step else 1
        mos = range(months.start, months.stop + step, step)

        if self.weatherconfig["file_granularity"] in {"daily", "dailymeans"}:
            # Separate files for each day (eg MERRA)
            # Check for complete set of files of year, month, day
            mo_tuples = [(yr, mo, monthrange(yr, mo)[1]) for yr in yrs for mo in mos]

            for mo_tuple in mo_tuples:
                # format: (yr, mo, number_days_in_month)
                yr, mo, nodays = mo_tuple
                for day in range(1, nodays + 1, 1):
                    filename = self.datasetfn(self.weatherconfig["fn"], yr, mo, day)
                    self.totalFiles.append((self.config, filename))
                    if not os.path.isfile(filename):
                        self.prepared = False
                        if check_complete:
                            logger.info("File `%s` not found!", filename)
                            incomplete_count += 1
                        if opendap and "url_opendap" in self.weatherconfig:
                            self.toDownload.append(
                                (
                                    self.config,
                                    filename,
                                    self.datasetfn_opendap(
                                        self.weatherconfig["url_opendap"],
                                        self.weatherconfig["variables"],
                                        yr,
                                        mo,
                                        day,
                                    ),
                                )
                            )
                        else:
                            if opendap:
                                logger.warning(
                                    "OpenDAP URL not specified for given dataset. "
                                    "Defaulting to standard download."
                                )

                            self.toDownload.append(
                                (
                                    self.config,
                                    filename,
                                    self.datasetfn(
                                        self.weatherconfig["url"], yr, mo, day
                                    ),
                                )
                            )
                    else:
                        self.downloadedFiles.append((self.config, filename))

        elif self.weatherconfig["file_granularity"] == "daily_multiple":
            # Separate files for each day (eg MERRA)
            # Check for complete set of files of year, month, day
            mo_tuples = [(yr, mo, monthrange(yr, mo)[1]) for yr in yrs for mo in mos]

            for mo_tuple in mo_tuples:
                # format: (yr, mo, number_days_in_month)
                yr, mo, nodays = mo_tuple
                for day in range(1, nodays + 1, 1):
                    filename = self.datasetfn(self.weatherconfig["fn"], yr, mo, day)
                    self.totalFiles.append((self.config, filename))
                    if not os.path.isfile(filename):
                        self.prepared = False
                        if check_complete:
                            logger.info("File `%s` not found!", filename)
                            incomplete_count += 1
                        if opendap and "url_opendap" in self.weatherconfig:
                            self.toDownload.append(
                                (
                                    self.config,
                                    filename,
                                    self.datasetfn_opendap(
                                        self.weatherconfig["url_opendap"][0],
                                        self.weatherconfig["variables_list"][0],
                                        yr,
                                        mo,
                                        day,
                                    ),
                                    self.datasetfn_opendap(
                                        self.weatherconfig["url_opendap"][1],
                                        self.weatherconfig["variables_list"][1],
                                        yr,
                                        mo,
                                        day,
                                    ),
                                )
                            )
                        else:
                            if opendap:
                                logger.warning(
                                    "OpenDAP URL not specified for given dataset. "
                                    "Defaulting to standard download."
                                )

                            self.toDownload.append(
                                (
                                    self.config,
                                    filename,
                                    self.datasetfn(
                                        self.weatherconfig["url"][0], yr, mo, day
                                    ),
                                    self.datasetfn(
                                        self.weatherconfig["url"][1], yr, mo, day
                                    ),
                                )
                            )
                    else:
                        self.downloadedFiles.append((self.config, filename))

        elif self.weatherconfig["file_granularity"] == "monthly":
            # Monthly files (eg ERA5)
            mo_tuples = [(yr, mo) for yr in yrs for mo in mos]
            for mo_tuple in mo_tuples:
                yr, mo = mo_tuple
                filename = self.datasetfn(self.weatherconfig["fn"], yr, mo)
                self.totalFiles.append((self.config, filename))
                if not os.path.isfile(filename):
                    self.prepared = False
                    if check_complete:
                        logger.info("File `%s` not found!", filename)
                        incomplete_count += 1
                    if self.module == "era5":
                        self.toDownload.append((self.config, filename, yr, mo))
                    else:
                        self.toDownload.append(
                            (
                                self.config,
                                filename,
                                self.datasetfn(self.weatherconfig["url"], yr, mo),
                            )
                        )
                else:
                    self.downloadedFiles.append((self.config, filename))

        elif self.weatherconfig["file_granularity"] == "monthly_multiple":
            mo_tuples = [(yr, mo) for yr in yrs for mo in mos]
            for mo_tuple in mo_tuples:
                yr, mo = mo_tuple
                filename = self.datasetfn(self.weatherconfig["fn"], yr, mo)
                self.totalFiles.append((self.config, filename))
                if not os.path.isfile(filename):
                    self.prepared = False
                    if check_complete:
                        logger.info("File `%s` not found!", filename)
                        incomplete_count += 1
                    self.toDownload.append(
                        (
                            self.config,
                            filename,
                            self.datasetfn(self.weatherconfig["url"][0], yr, mo),
                            self.datasetfn(self.weatherconfig["url"][1], yr, mo),
                        )
                    )
                else:
                    self.downloadedFiles.append((self.config, filename))

        if not self.prepared:
            if xs is None or ys is None:
                logger.warning(
                    "Arguments `xs` and `ys` not used in preparing dataset. Defaulting to global."
                )

            logger.info("%s files not completed.", incomplete_count)
        else:
            logger.info("Directory complete.")

    def datasetfn(self, fn, *args):
        """Construct file name from template function in `weather_data_config` and args (yr, mo, day)

        Args:
            fn (str): Template function in `weather_data_config`
            *args: Year, month, day

        Returns:
            str: File name
        """
        if len(args) == 3:
            dataset = dict(year=args[0], month=args[1], day=args[2])
        elif len(args) == 2:
            dataset = dict(year=args[0], month=args[1])
        else:
            return False
        if self.dataset_module.spinup_var:
            spinup = self.dataset_module.spinup_year(dataset["year"], dataset["month"])
            dataset.update({"spinup": spinup})

        return fn.format_map(dataset)

    def datasetfn_opendap(self, fn, variables, *args):
        """Construct url for OpenDap protocol. Depends on datasetfn()

        Args:
            fn (str): Template function in `weather_data_config`
            variables (list): List of variables to download
            *args: Year, month, day

        Returns:
            str: URL
        """

        # MERRA2 requires uppercase variable names
        if self.module == "merra2":
            variables = [v.upper() for v in variables]

        # TODO: Add lat lon subsetting based on self.bounds
        # opendap URL format requires knowing the indexes in MERRA2 -- not just lat, lon bounds --
        # e.g., ~/MERRA2/M2T1NXSLV.5.12.4/2021/01/MERRA2_400.tavg1_2d_slv_Nx.20210101.nc4.nc4?T2M[0:23][180:360][288:575],time,lat[180:360],lon[288:575]  # noqa: E501
        if self.bounds is not None:
            logger.info("Bounds not used in OpenDAP call. Defaulting to global.")

        fn = fn + "?" + ",".join(variables) + ",time,lat,lon"

        args = list(args)
        return self.datasetfn(fn, *args)

    def get_data(self, trim=True, testing=False):
        """Download data from server.

        Args:
            trim (bool, optional): Trim variables in file. Defaults to True.
            testing (bool, optional): Download only first file. Defaults to False.
        """

        if testing is True:
            if len(self.downloadedFiles) > 0:
                logger.info("First file for testing has already been downloaded.")
                return
            else:
                self.toDownload = [self.toDownload[0]]
                self.totalFiles = [self.totalFiles[0]]
                self.testDataset = True

        api_func = self.weatherconfig["api_func"]

        if self.module == "era5":
            api_func(
                self.toDownload,
                self.bounds,
                self.weatherconfig["keywords"],
                self.weatherconfig["product"],
                self.weatherconfig["product_type"],
                self.downloadedFiles,
            )

        elif self.module == "merra2":
            api_func(
                self.toDownload,
                self.weatherconfig["file_granularity"],
                self.downloadedFiles,
            )

        if self.downloadedFiles == self.totalFiles:
            self.prepared = True

        if trim is True and self.module not in config.untrimmable_datasets:
            self.trim_variables()

    def trim_variables(self):
        """Reduce size of file by trimming variables in file."""

        for f in self.downloadedFiles:
            file_path = f[1]
            variables = self.weatherconfig["variables"]

            with xr.open_dataset(file_path) as ds:
                var_rename = dict((v, v.lower()) for v in list(ds.data_vars))
                ds = ds.rename(var_rename)
                ds = ds[variables]

                fd, target = mkstemp(suffix=".nc4")
                os.close(fd)
                ds.to_netcdf(target)

            shutil.move(target, file_path)

    @property
    def meta_data_config(self):
        """Metadata configuration for dataset."""
        return dict(
            tasks_func=self.weatherconfig["tasks_func"],
            prepare_func=self.weatherconfig["meta_prepare_func"],
            template=self.weatherconfig["template"],
            file_granularity=self.weatherconfig["file_granularity"],
        )

    @property
    def weather_data_config(self) -> dict:
        """Weather data config for dataset."""
        return self.dataset_module.weather_data_config

    @property
    def projection(self) -> str:
        """Projection for dataset."""
        return self.dataset_module.projection

    # TODO: meta property/attribute does not exists here!
    # @property
    # def coords(self):
    #     """Coordinates for dataset."""
    #     return self.meta.coords  # pylint: disable=no-member

    @property
    def shape(self):
        """The shape of the Cutout by (y, x)."""
        return len(self.coords["y"]), len(self.coords["x"])

    @property
    def extent(self) -> Iterable[float]:
        """The extent of the Cutout by (x_min, x_max, y_min, y_max)."""
        return list(self.coords["x"].values[[0, -1]]) + list(
            self.coords["y"].values[[-1, 0]]
        )

    def grid_coordinates(self):
        """Return grid coordinates of the Cutout."""
        xs, ys = np.meshgrid(self.coords["x"], self.coords["y"])
        return np.asarray((np.ravel(xs), np.ravel(ys))).T

    def grid_cells(self):
        """Return grid cells of the Cutout."""
        coords = self.grid_coordinates()
        span = (coords[self.shape[1] + 1] - coords[0]) / 2
        return [box(*c) for c in np.hstack((coords - span, coords + span))]

    def __repr__(self):
        return "<Dataset {} years={}-{} months={}-{} datadir={} {}>".format(
            self.module,
            self.years.start,
            self.years.stop,
            self.months.start,
            self.months.stop,
            self.datadir,
            "prepared" if self.prepared else "unprepared",
        )
