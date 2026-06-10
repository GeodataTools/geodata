# Copyright 2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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
import os
import tempfile
from calendar import monthrange
from typing import Sequence

import numpy as np
import requests
import xarray as xr
from tqdm.auto import tqdm

from geodata.types import CoordRange, PathLike
from .._base import AtomicDataset, BaseDataset

logger = logging.getLogger(__name__)


def _convert_and_subset_lons_lats_merra2(
    ds: xr.Dataset, xs: CoordRange, ys: CoordRange
):
    if not isinstance(xs, slice):
        first, second, last = np.asarray(xs)[[0, 1, -1]]
        xs = slice(first - 0.1 * (second - first), last + 0.1 * (second - first))
    if not isinstance(ys, slice):
        first, second, last = np.asarray(ys)[[0, 1, -1]]
        ys = slice(first - 0.1 * (second - first), last + 0.1 * (second - first))

    ds = ds.sel(y=ys)

    # Lons should go from -180. to +180.
    if len(ds.coords["x"].sel(x=slice(xs.start + 360.0, xs.stop + 360.0))):
        ds = xr.concat(
            [ds.sel(x=slice(xs.start + 360.0, xs.stop + 360.0)), ds.sel(x=xs)], dim="x"
        )
        ds = ds.assign_coords(
            lon=np.where(
                ds.coords["x"].values <= 180,
                ds.coords["x"].values,
                ds.coords["x"].values - 360.0,
            )
        )
    else:
        ds = ds.sel(x=xs)

    return ds.assign_coords(lon=ds.coords["x"], lat=ds.coords["y"])


class MERRA2BaseDataset(BaseDataset):
    """MERRA2BaseDataset is a class that encaps a dataset from the MERRA2 reanalysis
    dataset. It provides a streamlined workflow for downloading, preprocessing,
    and storing of these datasets.
    """

    url_template: Sequence[str]
    variables: Sequence[str]
    module = "merra2"
    projection = "latlong"
    lat_direction = True
    frequency = "daily"

    def _download_file(self, file: AtomicDataset):
        assert hasattr(file, "url"), "URL is required to download the file"

        url: tuple[str] = file.url

        # Download the file
        match len(file.url):
            case 1:
                logger.debug("Downloading %s", url[0])
                with requests.get(url[0], stream=True) as r:
                    r.raise_for_status()

                    with (
                        open(file.path, "wb") as f,
                        tqdm(
                            total=int(r.headers.get("content-length", 0)),
                            unit="B",
                            unit_scale=True,
                        ) as pbar,
                    ):
                        for chunk in r.iter_content(chunk_size=8192):
                            f.write(chunk)
                            pbar.update(len(chunk))
            case _:
                with tempfile.TemporaryDirectory() as tempdir:
                    for i, url in enumerate(file.url):
                        logger.debug("Downloading %s", url)
                        with requests.get(url, stream=True) as r:
                            r.raise_for_status()

                            with (
                                open(os.path.join(tempdir, str(i)), "wb") as f,
                                tqdm(
                                    total=int(r.headers.get("content-length", 0)),
                                    unit="B",
                                    unit_scale=True,
                                ) as pbar,
                            ):
                                for chunk in r.iter_content(chunk_size=8192):
                                    f.write(chunk)
                                    pbar.update(len(chunk))

                    with xr.open_mfdataset(
                        [os.path.join(tempdir, str(i)) for i in range(len(file.url))],
                        combine="by_coords",
                        chunks="auto",
                    ) as ds:
                        ds.to_netcdf(file.path)
                        logger.info("Preprocessing complete with multiple files")

    def _dataset_postprocess(self, ds: xr.Dataset, **kwargs):
        """Postprocess the dataset after downloading and opening it.

        Args:
            ds (xr.Dataset): The dataset to postprocess.
            **kwargs: Additional keyword arguments.

        Returns:
            xr.Dataset: The postprocessed dataset.
        """

        # Rename the coordinates
        ds = ds.rename({"x": "lon", "y": "lat"})
        ds = ds.assign_coords(lon=ds.coords["lon"], lat=ds.coords["lat"])

        # Only keep the variables that are needed
        for var in ds.data_vars:
            if var.lower() not in self.variables:
                ds = ds.drop_vars(var)

        # Change all variables to lowercase
        return ds.rename({var: var.lower() for var in ds.data_vars})

    def spinup_year(self, year: int, month: int):
        """Returns the spinup period for the given year and month.
        See https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf for more
        information.

        Args:
            year (int): The year of the dataset
            month (int): The month of the dataset

        Returns:
            str: The spinup period
        """
        if year >= 1980 and year < 1992:
            spinup = "100"
        elif year >= 1992 and year < 2001:
            spinup = "200"
        elif year >= 2001 and year < 2011:
            spinup = "300"
        elif year >= 2011 and year < 2020:
            spinup = "400"
        elif year == 2020 and month == 9:
            spinup = "401"
        else:
            spinup = "400"

        return spinup

    def _daily_catalog(self):
        if not self.url_template:
            raise NotImplementedError("url_template is not defined for this dataset")

        catalog = super()._daily_catalog()

        for file in catalog:
            file.spinup = self.spinup_year(file.year, file.month)
            file.url = [t.format(**vars(file)) for t in self.url_template]

        return catalog

    def _monthly_catalog(self):
        if not self.url_template:
            raise NotImplementedError("url_template is not defined for this dataset")

        catalog = super()._monthly_catalog()

        for file in catalog:
            file.spinup = self.spinup_year(file.year, file.month)
            file.url = [t.format(**vars(file)) for t in self.url_template]

        return catalog

    @classmethod
    def meta_prepare_func(
        cls, xs: CoordRange, ys: CoordRange, year: int, month: int, **params
    ):
        with xr.open_mfdataset(cls._get_files(year, month), combine="by_coords") as ds:
            ds = ds.coords.to_dataset()
            ds = _convert_and_subset_lons_lats_merra2(ds, xs, ys)
            meta = ds.load()

        return meta

    @classmethod
    def tasks_func(
        cls,
        xs: CoordRange,
        ys: CoordRange,
        yearmonths: xr.DataArray,
        **meta_attrs,
    ):
        if not isinstance(xs, slice):
            xs = slice(*xs.values[[0, -1]])
        if not isinstance(ys, slice):
            ys = slice(*ys.values[[0, -1]])
        fn = meta_attrs["fn"]

        match cls.frequency:
            case "daily":
                logger.info(yearmonths)
                logger.info(
                    [
                        (year, month, day)
                        for year, month in yearmonths
                        for day in range(1, monthrange(year, month)[1] + 1, 1)
                    ]
                )

                return [
                    dict(
                        prepare_func=cls.prepare_func,
                        xs=xs,
                        ys=ys,
                        year=year,
                        month=month,
                        fn=fn.format(
                            year=year,
                            month=month,
                            day=day,
                            spinup=cls.spinup_year(year, month),
                        ),
                    )
                    for year, month in yearmonths
                    for day in range(1, monthrange(year, month)[1] + 1, 1)
                ]
            case "monthly":
                return [
                    dict(
                        xs=xs,
                        ys=ys,
                        year=year,
                        month=month,
                        fn=fn.format(year=year, month=month),
                    )
                    for year, month in yearmonths
                ]
            case _:
                raise NotImplementedError("Frequency not supported")

    @classmethod
    def prepare_func(
        cls,
        fn: PathLike,
        year: int,
        month: int,
        xs: slice,
        ys: slice,
        **kwargs,
    ):
        """Prepare the dataset for a given year and month."""

        if isinstance(fn, str) and not os.path.exists(fn):
            return
        if isinstance(fn, list) and not all(os.path.isfile(f) for f in fn):
            return

        with xr.open_dataset(fn) as ds:
            logger.info("Opening %s", fn)
            ds = _convert_and_subset_lons_lats_merra2(ds, xs, ys)
            ds = ds.rename({"x": "lon", "y": "lat"})
            ds = ds.assign_coords(lon=ds.coords["lon"], lat=ds.coords["lat"])

            yield (year, month), ds
