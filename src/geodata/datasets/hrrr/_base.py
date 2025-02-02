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
import pandas as pd
import xarray as xr

from .._base import BaseDataset

logger = logging.getLogger(__name__)


class HRRRBaseDataset(BaseDataset):
    """HRRRBaseDataset is a class that encaps a dataset from the HRRR
    dataset. It provides a streamlined workflow for downloading, preprocessing,
    and storing of these datasets.
    """

    module = "hrrr"
    projection = "latlong"
    frequency = "daily"
    _priority = ["google", "aws", "azure"]

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
        return ds.resample(time="1h").mean()
