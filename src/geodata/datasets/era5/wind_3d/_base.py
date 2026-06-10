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
import os

import xarray as xr

from geodata.types import PathLike
from .._base import ERA5BaseDataset, _subset_x_y_era5

logger = logging.getLogger(__name__)


class ERA5Wind3DBaseDataset(ERA5BaseDataset):
    """Base class for ERA5 3D wind datasets.
    
    This class provides the prepare_func implementation specific to wind_3d datasets,
    which use model levels and the reanalysis-era5-complete product.
    """

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
        
        This implementation is specific to wind_3d datasets which:
        - Use model levels (model_level coordinate)
        - Download from reanalysis-era5-complete product
        - Are stored as daily files
        """
        if isinstance(fn, str) and not os.path.exists(fn):
            return
        if isinstance(fn, list) and not all(os.path.isfile(f) for f in fn):
            return

        with xr.open_dataset(fn, engine="h5netcdf") as ds:
            logger.info("Opening %s", fn)
            ds = _subset_x_y_era5(ds, xs, ys)

            # New ERA5 format for hourly datasets
            # See https://forum.ecmwf.int/t/new-time-format-in-era5-netcdf-files/3796
            # TODO: We can remove this if we refactor geodata's convert module in the future
            if "valid_time" in ds.coords:
                ds = ds.rename({"valid_time": "time"})

            yield (year, month), ds

