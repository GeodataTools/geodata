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


import logging
import os.path as osp
import tempfile

import herbie
import numpy as np
import pandas as pd
import xarray as xr

from ...types import CoordRange
from .._base import BaseDataset

logger = logging.getLogger(__name__)


def _subset_x_y(ds: xr.Dataset, xs: slice, ys: slice):
    # Subset x,y according to xs, ys

    if not isinstance(xs, slice):
        first, second, last = np.asarray(xs)[[0, 1, -1]]
        xs = slice(first - 0.1 * (second - first), last + 0.1 * (second - first))
    if not isinstance(ys, slice):
        first, second, last = np.asarray(ys)[[0, 1, -1]]
        ys = slice(first - 0.1 * (second - first), last + 0.1 * (second - first))

    ds = ds.sel(y=ys)
    ds = ds.sel(x=xs)

    return ds


class HRRRBaseDataset(BaseDataset):
    """HRRRBaseDataset is a class that encaps a dataset from the HRRR
    dataset. It provides a streamlined workflow for downloading, preprocessing,
    and storing of these datasets.
    """

    module = "hrrr"
    projection = "latlong"
    frequency = "daily"
    _priority = ["google", "aws", "azure"]

    lat_direction = True

    def _extra_setup(self, **kwargs):
        self._herbie_save_dir = tempfile.TemporaryDirectory()

    def __del__(self):
        self._herbie_save_dir.cleanup()

    def _preprocess_individual_herbie(
        self, h: herbie.Herbie, search: str, hour: pd.DatetimeIndex
    ):
        tmp_ds = h.xarray(search, remove_grib=False)
        tmp_ds.to_netcdf(osp.join(self._herbie_save_dir.name, f"{hour}_{search}.nc"))
        tmp_ds.close()

        return osp.join(self._herbie_save_dir.name, f"{hour}_{search}.nc")

    def _dataset_postprocess(self, ds: xr.Dataset, **kwargs):
        # Because certain hours are missing, we need to reindex the dataset
        # to include all hours in the range

        logger.debug("Reindexing dataset to include all hours in the range")
        dt: np.datetime64 = (
            ds["time"].values[0].astype("datetime64[D]").astype("datetime64[m]")
        )
        ds = ds.reindex(time=pd.date_range(dt, periods=24, freq="h")).ffill(dim="time")

        # NOTE: For some reasons, HRRR's longitude is in the range of [0, 360]
        # instead of [-180, 180]. We need to put it back to [-180, 180].
        logger.debug("Fixing longitude range")
        ds["x"] = (ds["x"] % 360 + 540) % 360 - 180

        ds = ds.assign_coords(lon=ds.coords["x"], lat=ds.coords["y"])
        del ds["valid_time"]

        return ds

    @classmethod
    def meta_prepare_func(cls, xs: slice, ys: slice, year: int, month: int, **kwargs):
        with xr.open_mfdataset(cls._get_files(year, month), combine="by_coords") as ds:
            ds = ds.coords.to_dataset()
            ds = _subset_x_y(ds, xs, ys)
            meta = ds.load()

        return meta

    @classmethod
    def tasks_func(
        cls,
        xs: CoordRange,
        ys: CoordRange,
        yearmonths: xr.DataArray,
        **kwargs,
    ):
        if not isinstance(xs, slice):
            xs = slice(*xs.values[[0, -1]])
        if not isinstance(ys, slice):
            ys = slice(*ys.values[[0, -1]])

        return [
            dict(
                prepare_func=cls.prepare_func,
                xs=xs,
                ys=ys,
                year=year,
                month=month,
                fn=cls._get_files(year, month, **kwargs),
            )
            for year, month in yearmonths
        ]

    @staticmethod
    def prepare_func(fn, year, month, xs, ys, **kwargs):
        if isinstance(fn, str) and not osp.exists(fn):
            return
        if isinstance(fn, list) and not all(osp.isfile(f) for f in fn):
            return

        with xr.open_mfdataset(fn) as ds:
            logger.info("Opening %s", fn)

            yield (year, month), _subset_x_y(ds, xs, ys)
