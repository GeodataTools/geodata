# Copyright 2024-2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD), Keyu Long (UCSD)

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

import numpy as np
import xarray as xr

from geodata.types import CoordRange, PathLike
from .._base import BaseDataset

logger = logging.getLogger(__name__)

try:
    import cdsapi
except ImportError:
    logger.warning(
        "cdsapi is not installed. You will not be able to download ERA5 data."
        "Please install it with `pip install 'geodata-re[download]'."
    )
    cdsapi = None


def _convert_and_subset_lons_lats_era5(ds: xr.Dataset, xs: slice, ys: slice):
    # Longitudes should go from -180. to +180.
    if len(ds.coords["x"].sel(x=slice(xs.start + 360.0, xs.stop + 360.0))):
        ds = xr.concat(
            [ds.sel(x=slice(xs.start + 360.0, xs.stop + 360.0)), ds.sel(x=xs)], dim="x"
        )
        ds = ds.assign_coords(
            x=np.where(
                ds.coords["x"].values <= 180,
                ds.coords["x"].values,
                ds.coords["x"].values - 360.0,
            )
        )

    # Subset x and y
    return _subset_x_y_era5(ds, xs, ys)


def _subset_x_y_era5(ds: xr.Dataset, xs: slice, ys: slice):
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


class ERA5BaseDataset(BaseDataset):
    """ERA5BaseDataset is a class that encaps a dataset from the ERA5 reanalysis
    dataset. It provides a streamlined workflow for downloading, preprocessing,
    and storing of these datasets.
    """

    module = "era5"
    projection = "latlong"
    lat_direction = False

    def _extra_setup(self, **kwargs):
        self.logger = logging.getLogger(__name__.replace("._base", ".client"))
        self.client = cdsapi.Client(
            info_callback=self.logger.info,
            error_callback=self.logger.error,
            debug_callback=self.logger.debug,
            warning_callback=self.logger.warning,
        )

    @classmethod
    def meta_prepare_func(cls, xs: slice, ys: slice, year: int, month: int, **kwargs):
        # Reference of the quantities
        # https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
        # Geopotential is aka Orography in the CDS:
        # https://confluence.ecmwf.int/pages/viewpage.action?pageId=78296105

        with xr.open_mfdataset(cls._get_path(year, month), combine="by_coords") as ds:
            ds = ds.coords.to_dataset()
            ds = _convert_and_subset_lons_lats_era5(ds, xs, ys)
            meta = ds.load()

        return meta

    @classmethod
    def tasks_func(
        cls,
        xs: CoordRange,
        ys: CoordRange,
        yearmonths: xr.DataArray,
        prepare_func: callable,
        **meta_attrs,
    ):
        if not isinstance(xs, slice):
            xs = slice(*xs.values[[0, -1]])
        if not isinstance(ys, slice):
            ys = slice(*ys.values[[0, -1]])
        fn = meta_attrs["fn"]

        logger.info(yearmonths)
        logger.info(list(yearmonths))

        return [
            dict(
                prepare_func=prepare_func,
                xs=xs,
                ys=ys,
                year=year,
                month=month,
                fn=fn.format(year=year, month=month),
            )
            for year, month in yearmonths
        ]

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
        """Prepare the dataset for a given year and month.
        
        This method should be overridden by subclasses (e.g., ERA5Wind3DBaseDataset,
        ERA5WindSolarBaseDataset) to provide dataset-specific preparation logic.
        """
        raise NotImplementedError(
            "prepare_func must be implemented by a subclass. "
            "Use ERA5Wind3DBaseDataset or ERA5WindSolarBaseDataset, "
            "or override this method in your subclass."
        )

    def _dataset_postprocess(self, ds, **kwargs):
        return super()._dataset_postprocess(ds, **kwargs)
